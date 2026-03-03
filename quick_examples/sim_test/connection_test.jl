# =============================================================================
# connection_test.jl
# Trollheim HPP — Component-Based (Acausal) Connection Model
# GGOV1 Droop Governor, 10 % load step at t = 5 s
#
# This file demonstrates the full acausal component architecture:
#   Reservoir → Penstock → FrancisTurbineAffinity → Reservoir (hydraulic)
#   FrancisTurbineAffinity → RotorInertia → SimpleGenerator (mechanical)
#   RotorInertia + FrancisTurbineAffinity → GGOV1Governor → FrancisTurbineAffinity (control)
#
# Connection strategy
#   Physical domains (HydraulicPort, RotationalPort) : connect()
#   Signal domains (speed, power, gate)              : algebraic equations
#     governor.speed_in.u ~ rotor.omega
#     governor.power_in.u ~ turbine.P_mech
#     turbine.opening.u   ~ governor.gate_out.u
#
# Note: RotorInertia and FrancisTurbineAffinity intentionally have no SignalOutPorts
#       (acausal physical ports only). Signal extraction is done algebraically
#       in the composite @equations block — a standard MTK pattern.
#
# Operating Point B:  H = 371 m,  Q = 37 m³/s,  P ≈ 130 MW,  n = 375 RPM
#
# Output:  images/connection_test.png
# =============================================================================

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots
using Printf

# ── Operating Point B ─────────────────────────────────────────────────────────
const H_n       = 371.0
const Q_n       = 37.0
const P_r       = 130.0e6             # rated power [W]
const n_rpm     = 375.0
const omega_ref = 2π * n_rpm / 60    # ≈ 39.27 rad/s
const ρ         = 1000.0
const g_acc     = 9.81

# ── Penstock geometry (Trollheim benchmark, two-section penstock) ─────────────
#   Section 1: L=363 m, D=4.7 m    Section 2: L=145 m, D=3.3 m
#   Lump into single equivalent penstock: L=508 m, D=4.1 m (area-weighted)
const L_p   = 363.0 + 145.0          # total penstock length [m]
const D_p   = 4.1                    # equivalent diameter [m]

# ── Turbine sizing (OP-B) ─────────────────────────────────────────────────────
const D_runner = 2.5
const τ_o_ss   = 0.90
const K_q_v    = Q_n / (τ_o_ss * D_runner^2 * sqrt(H_n))   # ≈ 0.3415
const η_max_v  = 0.97
const c_eta_v  = 0.25
const Q_rated_v = Q_n

# ── Rotor: mechanical starting time M_a = 8 s ────────────────────────────────
const J_rotor  = 8.0 * P_r / omega_ref^2   # ≈ 6.74×10⁵ kg·m²

# ── Governor: GGOV1 droop (Trollheim benchmark σ = 0.10) ─────────────────────
const sigma   = 0.10      # static droop (dimensionless)
const T_ps    = 1.75      # governor lag [s]  (= T_gov)
const T_gs    = 0.20      # gate servo time constant [s]
const K_i_gov = 1.0 / T_gs   # gate integrator gain = 5.0

# GGOV1 uses:  P_ref = P_set + (omega_ref - omega) / R_droop
# For dimensional consistency (P_set in W, omega in rad/s):
#   R_droop = sigma * omega_ref / P_r   [rad·s·W⁻¹]
const R_droop_phy = sigma * omega_ref / P_r   # ≈ 3.02e-8 [rad/s/W]

println("=== Trollheim HPP — Component-Based GGOV1 Droop Simulation ===")
@printf "  ω_ref     = %.4f rad/s  (%g RPM)\n" omega_ref n_rpm
@printf "  K_q       = %.4f\n"                  K_q_v
@printf "  J_rotor   = %.3e kg·m²\n"            J_rotor
@printf "  sigma     = %.2f   T_gov = %.2f s   K_i = %.2f\n" sigma T_ps K_i_gov
@printf "  R_droop   = %.4e  (= sigma * omega_ref / P_r)\n"   R_droop_phy
println()

# =============================================================================
# ── Composite component model ─────────────────────────────────────────────────
# =============================================================================
@mtkmodel TrollheimDroop begin
    @components begin
        # ── Hydraulic circuit ─────────────────────────────────────────────────
        reservoir = Reservoir(
            H     = H_n,
            rho   = ρ,
            g     = g_acc,
        )
        penstock = Penstock(
            L      = L_p,
            D_pipe = D_p,
            rho = ρ,
            e   = 1.5e-5,           # commercial steel roughness [m]
        )
        turbine = FrancisTurbineAffinity(
            D       = D_runner,
            K_q     = K_q_v,
            eta_max = η_max_v,
            c_eta   = c_eta_v,
            Q_rated = Q_rated_v,
            rho     = ρ,
            g       = g_acc,
        )
        tailwater = Reservoir(
            H     = 0.0,
            rho   = ρ,
            g     = g_acc,
        )

        # ── Mechanical train ──────────────────────────────────────────────────
        rotor = RotorInertia(
            J       = J_rotor,
            b_v     = 50.0,         # viscous bearing damping [N·m·s/rad]
            omega_0 = omega_ref,
        )
        generator = SimpleGenerator(
            P_rated     = P_r,
            omega_rated = omega_ref,
            omega_s     = omega_ref,
            D_d         = 1.5,      # synchronising damping
            eta_gen     = 0.985,
        )

        # ── GGOV1 droop governor ──────────────────────────────────────────────
        governor = GGOV1Governor(
            R_droop   = R_droop_phy,   # sigma * omega_ref / P_r → physical units
            T_p       = 0.05,
            T_gov     = T_ps,
            omega_ref = omega_ref,
            P_set     = P_r,        # rated load setpoint [W]
            tau_min   = 0.0,
            tau_max   = 1.0,
            K_i       = K_i_gov,
            tau_0     = τ_o_ss,
        )
    end

    @equations begin
        # ── Hydraulic: connect() on HydraulicPort ─────────────────────────────
        connect(reservoir.port,    penstock.port_a)
        connect(penstock.port_b,   turbine.port_in)
        connect(turbine.port_out,  tailwater.port)

        # ── Mechanical: connect() on RotationalPort ───────────────────────────
        connect(turbine.shaft,          rotor.flange_turbine)
        connect(rotor.flange_generator, generator.flange)

        # ── Control signals: algebraic equations (no SignalOutPort on rotor/turbine)
        # Speed feedback:  governor reads rotor angular velocity
        governor.speed_in.u ~ rotor.omega

        # Power feedback:  governor reads turbine mechanical power
        governor.power_in.u ~ turbine.P_mech

        # Gate command:    governor drives turbine guide vane opening
        turbine.opening.u   ~ governor.gate_out.u
    end
end

# ── Compile (MTK structural simplification) ───────────────────────────────────
@mtkcompile sys = TrollheimDroop()

# ── Initial conditions ────────────────────────────────────────────────────────
# Steady state at gate = 0.90, rated speed, P_meas initialised at SS power
Q_ic    = τ_o_ss * K_q_v * D_runner^2 * sqrt(H_n)       # ≈ 37 m³/s
P_ic    = ρ * g_acc * Q_ic * H_n * η_max_v               # ≈ 130 MW
dm_ic   = ρ * Q_ic                                        # ≈ 37 000 kg/s
η_ic    = η_max_v * (1.0 - c_eta_v * (Q_ic / Q_rated_v - 1.0)^2)

# Penstock algebraic ICs
A_pipe  = π * D_p^2 / 4
μ_water = 1e-3                              # dynamic viscosity [Pa·s]
Re_ic   = dm_ic * D_p / (μ_water * A_pipe) # ≈ 4×10⁷ (fully turbulent)
f_D_ic  = 1.0 / (2.0 * log10(1.5e-5 / (3.7 * D_p) + 5.7 / Re_ic^0.9))^2  # Swamee-Jain
dp_f_ic = f_D_ic * (L_p / D_p) * (dm_ic / (ρ * A_pipe))^2 * ρ / 2

u0 = Dict(
    sys.rotor.omega          => omega_ref,
    sys.rotor.theta          => 0.0,
    sys.penstock.dm          => dm_ic,
    sys.penstock.p_avg       => ρ * g_acc * H_n,
    # MTSL block states: pow_lag / gov_lag (FirstOrder) and gate_int (Integrator)
    # P_meas / x_gov / x_int are now algebraic alias variables — set the .x states.
    sys.governor.pow_lag.x   => P_ic,    # pow_lag output = k*x = P_meas at t=0
    sys.governor.gov_lag.x   => 0.0,     # gov_lag output = x_gov at t=0
    sys.governor.gate_int.x  => τ_o_ss,  # integrator state = x_int at t=0
)

# ── Algebraic guesses — break MTK initialization substitution cycles ──────────
# All purely algebraic (non-differential) variables need a numeric seed so MTK
# can solve the DAE initialization without entering a symbolic cycle.
guesses = Dict(
    # Turbine algebraics
    sys.turbine.H         => H_n,
    sys.turbine.Q         => Q_ic,
    sys.turbine.eta       => η_ic,
    sys.turbine.P_mech    => P_ic,
    sys.turbine.dm        => dm_ic,
    sys.turbine.tau_shaft => P_ic / omega_ref,
    # Penstock algebraics
    sys.penstock.Re       => Re_ic,
    sys.penstock.f_D      => f_D_ic,
    sys.penstock.dp_f     => dp_f_ic,
)

tspan = (0.0, 60.0)
prob  = ODEProblem(sys, u0, tspan; guesses = guesses)

# ── 10 % load step at t = 5 s — step governor P_set (power demand setpoint) ──
cb = DiscreteCallback(
    (u, t, i) -> t == 5.0,
    i -> (i.ps[sys.governor.P_set] = 1.1 * P_r; nothing),
    save_positions = (true, true)
)

# ── Solve ─────────────────────────────────────────────────────────────────────
println("Solving...")
@time sol = solve(prob, Rodas5P(),
                  tstops    = [5.0],
                  callback  = cb,
                  abstol    = 1e-7,
                  reltol    = 1e-7)

println("  Retcode: ", sol.retcode)
ω_ss = sol(60.0; idxs = sys.rotor.omega)
Δω_ss = ω_ss - omega_ref
# Theoretical SS deviation: Δω_ss = (ΔP / P_r) × sigma × omega_ref
# For 10% step: ΔP/P_r = 0.1 → Δω_ss = -0.1 × sigma × omega_ref
Δω_pred = -0.1 * sigma * omega_ref    # theoretical droop drop for 10% step
@printf "  ω_ss (t=60)     = %.5f rad/s   (ref = %.5f)\n" ω_ss omega_ref
@printf "  Δω_ss measured  = %+.5f rad/s\n" Δω_ss
@printf "  Δω_ss predicted = %+.5f rad/s  (σ=%.2f, 10%% load step)\n" Δω_pred sigma
println()

# ── Plot ──────────────────────────────────────────────────────────────────────
default(linewidth = 1.8, tickfontsize = 9, guidefontsize = 10,
        legendfontsize = 8, left_margin = 8Plots.mm, bottom_margin = 3Plots.mm)

tv = range(0.0, 60.0, 600)

ω_v    = [sol(t; idxs = sys.rotor.omega)           for t in tv]
τo_v   = [sol(t; idxs = sys.governor.tau_o)        for t in tv]
Pm_v   = [sol(t; idxs = sys.turbine.P_mech) / 1e6  for t in tv]
Δω_v   = [(ω - omega_ref) / omega_ref * 100        for ω in ω_v]

p_omega = plot(tv, ω_v;
               label  = "ω [rad/s]",  color  = :steelblue,
               ylabel = "Speed  [rad/s]",
               title  = "Trollheim HPP — Component Model  |  GGOV1 Droop (σ=0.10)  |  10% step at t=5 s")
hline!([omega_ref];                ls = :dash, lc = :black,  lw = 1, label = "ω_ref")
hline!([omega_ref + Δω_pred];      ls = :dot,  lc = :orange, lw = 1, label = "Δω_ss predicted")
vline!([5.0];                      ls = :dot,  lc = :grey,   lw = 1, label = "step")

p_gate = plot(tv, τo_v;
              label  = "τ_o [-]", color = :steelblue,
              ylabel = "Gate opening  [pu]", legend = :bottomright)
vline!([5.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p_power = plot(tv, Pm_v;
               label  = "P_mech [MW]", color = :steelblue,
               ylabel = "Mechanical power  [MW]", legend = :topleft)
hline!([P_r / 1e6];       ls = :dash, lc = :black,  lw = 1, label = "P_rated")
hline!([1.1 * P_r / 1e6]; ls = :dot,  lc = :orange, lw = 1, label = "1.1 P_rated (step target)")
vline!([5.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p_freq = plot(tv, Δω_v;
              label  = "Δω/ω_ref [%]", color = :steelblue,
              ylabel = "Freq. deviation  [%]", legend = :bottomleft)
hline!([Δω_pred / omega_ref * 100]; ls = :dot, lc = :orange, lw = 1,
       label = "Predicted Δω_ss (droop)")
hline!([0.0];                       ls = :dash, lc = :black, lw = 1, label = "0 %")
vline!([5.0]; ls = :dot, lc = :grey, lw = 1, label = false)

ptop = plot(p_omega, p_gate, p_power, p_freq;
            layout        = (4, 1),
            size          = (950, 1100),
            link          = :x,
            xlabel        = "Time  [s]")

mkpath(joinpath(@__DIR__, "images"))
out_png = joinpath(@__DIR__, "images", "connection_test.png")
savefig(ptop, out_png)
println("Plot saved  ->  sim_test/images/connection_test.png")
