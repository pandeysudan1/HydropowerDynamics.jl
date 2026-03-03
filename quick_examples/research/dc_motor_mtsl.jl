# =============================================================================
# dc_motor_mtsl.jl
# Research: ModelingToolkitStandardLibrary.jl DC Motor → Hydropower Analogy
#
# PURPOSE
# -------
# Study how MTSL handles the three connection kinds that also appear in
# connection_test.jl, and map each pattern to the hydropower equivalents.
#
# THREE CONNECTION KINDS IN MTSL
# ──────────────────────────────
#  1. ACAUSAL physical ports (Kirchhoff-type)
#       Electrical  : voltage (across) + current (through)
#       Rotational  : angular velocity (across) + torque (through)
#       Hydraulic   : pressure p (across) + mass-flow dm (through)
#     → connected with  connect(a.port, b.port)
#     → MTK generates  Kirchhoff sum equations automatically
#
#  2. CAUSAL "signal" ports  (one-directional data flow)
#       RealInput / RealOutput   –  just a single variable u
#     → connected with  connect(source.output, sink.input)
#       or named analysis points: connect(sensor.w, :y, controller.input2)
#
#  3. CROSS-DOMAIN bridges  (physical ↔ signal)
#       SpeedSensor  : RotationalPort flange  →  SignalOutPort  w
#       EMF          : ElectricalPort + RotationalPort (coupled: V = k·ω, τ = k·i)
#     → these are the only components that cross the acausal/causal boundary
#
# HYDROPOWER ANALOGY TABLE
# ────────────────────────
#  DC Motor domain              │ Hydropower domain
#  ─────────────────────────────┼────────────────────────────────────────────
#  Voltage source (Voltage)     │ Reservoir   (fixed-pressure hydraulic port)
#  Resistor R1                  │ Penstock friction loss  (dp ~ dm²)
#  Inductor L1                  │ Penstock water inertia  (D(dm) ~ dp/L)
#  EMF                          │ FrancisTurbine  (hydraulic ↔ mechanical)
#  Inertia J                    │ RotorInertia
#  Damper f                     │ Bearing damping (b_v)
#  SpeedSensor                  │ (missing in HydroPowerDynamics!) → need to add
#  Feedback (error = ref − meas)│ ISO/GGOV1 error signal
#  LimPI  (PI + anti-windup)   │ ISO PID governor (Kp, Ki, rate limits)
#
# KEY DIFFERENCE FROM connection_test.jl
# ─────────────────────────────────────────
# In connection_test.jl, the speed signal is extracted ALGEBRAICALLY:
#     governor.speed_in.u ~ rotor.omega
# The MTSL pattern uses a proper SpeedSensor component which has a
# RotationalPort flange (acausal) and outputs a causal signal:
#     connect(inertia.flange_b, speed_sensor.flange)
#     connect(speed_sensor.w, :y, feedback.input2)
# The ":y" is an "analysis point" — a named tap that also enables
# automatic sensitivity/loop-transfer-function computation.
#
# WHY ANALYSIS POINTS MATTER FOR HYDROPOWER CONTROL
# ──────────────────────────────────────────────────
# Once declared, analysis points allow:
#     Blocks.get_sensitivity(model, :y, op = ...)   → S(s) = 1/(1+PC)
#     Blocks.get_comp_sensitivity(model, :y, ...)   → T(s) = PC/(1+PC)
#     Blocks.get_looptransfer(model, :y, ...)       → L(s) = PC
# This is invaluable for ISO/GGOV1 droop tuning: verify GM, PM, bandwidth.
#
# Output:  quick_examples/research/images/dc_motor_mtsl.png
# =============================================================================

using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using Plots
using Printf

# =============================================================================
# ── 1. DC MOTOR MODEL (direct from MTSL tutorial) ─────────────────────────────
#
# Circuit:  Voltage → R1 → L1 → EMF → Ground
# Mechanical: EMF flange ─ friction ─ inertia ─ load torque
# Control:  SpeedSensor → Feedback → LimPI → Voltage.V
#
# Analysis points are tagged with :y (plant output) and :u (plant input)
# so sensitivity functions can be extracted directly.
# =============================================================================
@mtkmodel DCMotor begin
    @parameters begin
        R         = 0.5,   [description = "Armature resistance [Ω]"]
        L         = 4.5e-3,[description = "Armature inductance [H]"]
        k         = 0.5,   [description = "Motor constant [N·m/A = V·s/rad]"]
        J         = 0.02,  [description = "Rotor inertia [kg·m²]"]
        f         = 0.01,  [description = "Viscous friction [N·m·s/rad]"]
        tau_L_step = -0.3, [description = "Load torque step amplitude [N·m]"]
    end
    @components begin
        # ── Electrical ────────────────────────────────────────────────────────
        ground        = Ground()
        source        = Voltage()                            # voltage actuator
        R1            = Resistor(R = R)
        L1            = Inductor(L = L)
        emf           = EMF(k = k)                          # electro-mechanical bridge

        # ── Mechanical ───────────────────────────────────────────────────────
        fixed         = Fixed()                             # ground for rotation
        friction      = Damper(d = f)
        inertia       = Inertia(J = J)
        load          = Torque()
        speed_sensor  = SpeedSensor()                       # acausal → causal bridge

        # ── Control signals (causal, Blocks domain) ──────────────────────────
        ref           = Blocks.Step(height = 1, start_time = 0)   # speed setpoint [rad/s]
        load_step     = Blocks.Step(height = tau_L_step, start_time = 3)
        feedback      = Blocks.Feedback()                  # error = ref − measured
        pi_controller = Blocks.LimPI(k = 1.1, T = 0.035, u_max = 10, Ta = 0.035)
        # LimPI: PI with anti-windup limiter
        #   k    = proportional gain
        #   T    = integral time constant  (Ki = k/T)
        #   u_max= output saturation       (anti-windup)
        #   Ta   = anti-windup time constant
    end
    @equations begin
        # ── Electrical circuit (ACAUSAL) ──────────────────────────────────────
        # MTK generates KVL automatically: sum of across (voltages) = 0
        # and KCL: sum of through (currents) = 0 at each node.
        connect(source.p, R1.p)
        connect(R1.n,     L1.p)
        connect(L1.n,     emf.p)
        connect(emf.n,    source.n, ground.g)

        # ── Mechanical (ACAUSAL) ──────────────────────────────────────────────
        # MTK generates torque-balance automatically: sum of torques = 0
        # and angular velocity equality at connected flanges.
        connect(fixed.flange,          emf.support, friction.flange_b)
        connect(emf.flange,            friction.flange_a, inertia.flange_a)
        connect(inertia.flange_b,      load.flange)
        connect(inertia.flange_b,      speed_sensor.flange)

        # ── Cross-domain bridges ──────────────────────────────────────────────
        # SpeedSensor: RotationalPort → SignalOutPort  w  [rad/s]
        # load step:   causal Blocks.Step → Torque component
        connect(load_step.output,      load.tau)

        # ── Control loop (CAUSAL signals) ────────────────────────────────────
        # Analysis point :y  at plant output (measured speed)
        # Analysis point :u  at plant input  (voltage command)
        # These let get_sensitivity(), get_looptransfer() extract L(s),S(s),T(s)
        connect(ref.output,            feedback.input1)
        connect(speed_sensor.w, :y,    feedback.input2)     # tagged → analysis
        connect(feedback.output,       pi_controller.err_input)
        connect(pi_controller.ctr_output, :u, source.V)     # tagged → analysis
    end
end

# ── Compile ───────────────────────────────────────────────────────────────────
@named model = DCMotor()
sys = mtkcompile(model)

println("=== DC Motor — MTSL acausal+causal simulation ===")
println()
println("Structural info:")
println("  States      : ", length(unknowns(sys)))
println("  Parameters  : ", length(ModelingToolkit.parameters(sys)))
println()

# ── Solve ─────────────────────────────────────────────────────────────────────
# Only L1.i needs an explicit IC (the others are determined by initialization).
prob = ODEProblem(sys, [sys.L1.i => 0.0], (0.0, 6.0))
println("Solving (0 → 6 s)...")
@time sol = solve(prob, Rodas5P())
println("  Retcode: ", sol.retcode)
println()

# ── Print steady-state ────────────────────────────────────────────────────────
ω_ss  = sol(1.9; idxs = sys.inertia.w)   # just before load step
ω_end = sol(5.9; idxs = sys.inertia.w)   # after load step settled
@printf "  ω_ss  (t=1.9 s, pre-step)  = %.4f rad/s   (setpoint = 1.0)\n" ω_ss
@printf "  ω_end (t=5.9 s, post-step) = %.4f rad/s\n"                   ω_end
println()

# =============================================================================
# ── 2. PLOT ────────────────────────────────────────────────────────────────────
# =============================================================================
default(linewidth = 1.8, tickfontsize = 9, guidefontsize = 10,
        legendfontsize = 8, left_margin = 8Plots.mm, bottom_margin = 3Plots.mm)

tv = range(0.0, 6.0, 600)

ω_v     = [sol(τ; idxs = sys.inertia.w)     for τ in tv]
i_v     = [sol(τ; idxs = sys.L1.i)          for τ in tv]
V_v     = [sol(τ; idxs = sys.source.V.u)    for τ in tv]
τ_L_v   = [sol(τ; idxs = sys.load.tau.u)    for τ in tv]
ref_v   = [sol(τ; idxs = sys.ref.output.u)  for τ in tv]

p1 = plot(tv, ω_v;   label = "ω meas [rad/s]", color = :steelblue,
          ylabel = "Speed  [rad/s]",
          title  = "DC Motor — MTSL PI Speed Control  |  load step at t=3 s")
plot!(p1, tv, ref_v; label = "ω ref",    color = :black,   ls = :dash, lw = 1)
vline!([3.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p2 = plot(tv, i_v;   label = "i_a [A]",         color = :steelblue,
          ylabel = "Armature current  [A]", legend = :topright)
vline!([3.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p3 = plot(tv, V_v;   label = "V_source [V]",     color = :steelblue,
          ylabel = "Input voltage  [V]", legend = :topright)
vline!([3.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p4 = plot(tv, τ_L_v; label = "τ_load [N·m]",    color = :firebrick,
          ylabel = "Load torque  [N·m]", legend = :bottomleft)
vline!([3.0]; ls = :dot, lc = :grey, lw = 1, label = false)

ptop = plot(p1, p2, p3, p4;
            layout        = (4, 1),
            size          = (950, 1100),
            link          = :x,
            xlabel        = "Time  [s]")

mkpath(joinpath(@__DIR__, "images"))
out_png = joinpath(@__DIR__, "images", "dc_motor_mtsl.png")
savefig(ptop, out_png)
println("Plot saved  →  research/images/dc_motor_mtsl.png")
println()

# =============================================================================
# ── 3. MAPPING SUMMARY  (printed to terminal)
# =============================================================================
println("""
╔═══════════════════════════════════════════════════════════════════════════╗
║  PATTERN ANALYSIS: DC Motor  ↔  HydropowerDynamics  connection_test.jl  ║
╠═══════════════════════════════════════════════════════════════════════════╣
║                                                                           ║
║  A. ACAUSAL PHYSICAL CONNECTIONS  (both projects use connect())           ║
║  ─────────────────────────────────────────────────────────────────────    ║
║  DC Motor              │ HydroPowerDynamics                               ║
║  ─────────────────────────────────────────────────────────────────────    ║
║  Electrical domain:    │ Hydraulic domain:                                ║
║    across = voltage    │   across = pressure  p  [Pa]                    ║
║    through = current   │   through = mass-flow dm [kg/s]                 ║
║  Rotational domain:    │ Rotational domain (same!):                       ║
║    across = w [rad/s]  │   across = omega  [rad/s]                       ║
║    through = tau [N·m] │   through = tau   [N·m]                         ║
║                        │                                                  ║
║  B. CAUSAL SIGNAL CONNECTIONS                                             ║
║  ─────────────────────────────────────────────────────────────────────    ║
║  Blocks.Step            → speed reference                                 ║
║  Blocks.Feedback        → error  e = ω_ref - ω                           ║
║  Blocks.LimPI           → PI controller  (replaces custom ISO PID)       ║
║  → governor.speed_in.u ~ rotor.omega is the ALGEBRAIC equivalent         ║
║    but MTSL's SpeedSensor + connect() is the proper acausal way          ║
║                                                                           ║
║  C. CROSS-DOMAIN BRIDGES                                                  ║
║  ─────────────────────────────────────────────────────────────────────    ║
║  DC Motor: EMF  (V=k·ω, τ=k·i)  →  bridges electrical ↔ mechanical      ║
║  Hydro:    FrancisTurbineAffinity   →  bridges hydraulic ↔ mechanical    ║
║                                                                           ║
║  DC Motor: SpeedSensor (flange → signal w)  →  acausal→causal            ║
║  Hydro:    currently MISSING – using algebraic hack                       ║
║            IMPROVEMENT: add RotationalSpeedSensor to src/mechanical.jl   ║
║            then: connect(rotor.flange_gen, speed_sensor.flange)          ║
║                  connect(speed_sensor.w, :y, governor.speed_in)          ║
║                                                                           ║
║  D. ANALYSIS POINTS  (unique MTSL feature)                               ║
║  ─────────────────────────────────────────────────────────────────────    ║
║  Named tap in connect():  connect(sensor.w, :y, controller.input2)       ║
║  Enables:                                                                 ║
║    get_sensitivity(model, :y, op=...)   → S(s) = 1/(1+PC)               ║
║    get_comp_sensitivity(model, :y, ...) → T(s)                           ║
║    get_looptransfer(model, :y, ...)     → L(s) = PC                      ║
║  → ControlSystemsBase.bodeplot(So, To)  shows gain/phase margin          ║
║  HYDROPOWER APPLICATION: tune droop σ and ISO Kp/Ki with bode plots     ║
║                                                                           ║
║  E. HYDRAULIC LIBRARY IN MTSL                                             ║
║  ─────────────────────────────────────────────────────────────────────    ║
║  MTSL has ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible║
║    Tube         – penstock with inertia + built-in friction_factor        ║
║    FixedPressure– reservoir (fixed pressure boundary)                    ║
║    Valve        – guide vane with area input signal                       ║
║    Volume       – surge tank (moving-wall hydraulic ↔ mechanical)        ║
║  → These could eventually replace custom hydraulic.jl components        ║
║    but they use gauge pressure and isothermal density, needing            ║
║    adaptation for water-hammer (large B) and surge tank gravity head.    ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
""")
