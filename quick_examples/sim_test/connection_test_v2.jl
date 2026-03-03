# =============================================================================
# connection_test_v2.jl
# Trollheim HPP — Component-Based Model  (v2: sensors + analysis points)
#
# WHAT IS DIFFERENT FROM connection_test.jl (v1)
# ───────────────────────────────────────────────
#  v1  signal wiring used algebraic equation hacks:
#       governor.speed_in.u ~ rotor.omega          ← no analysis point possible
#       turbine.opening.u   ~ governor.gate_out.u  ← no analysis point possible
#
#  v2  replaces those with proper sensor components + named analysis points:
#       connect(rotor.flange_generator, speed_sensor.flange, generator.flange)
#       connect(speed_sensor.w,         :y_omega,  governor.speed_in)
#       connect(governor.gate_out,      :u_gate,   turbine.opening)
#
#  The :y_omega / :u_gate tags make the closed-loop linearisable:
#       Blocks.get_sensitivity(model, :y_omega; op)  →  S(s) = 1/(1+PC)
#       Blocks.get_looptransfer(model, :y_omega; op) →  L(s) = PC
#       ControlSystemsBase.margin(L)                 →  GM, PM
#
# ARCHITECTURE
#   Hydraulic  : Reservoir → Penstock → FrancisTurbineAffinity → Reservoir
#   Mechanical : turbine.shaft → RotorInertia → SimpleGenerator
#   Sensors    : RotationalSpeedSensor (parallel tap on flange_generator node)
#   Governor   : GGOV1Governor (droop, power lag, gate integrator)
#   Loop       : speed_sensor.w :y_omega→ governor.speed_in
#                governor.gate_out :u_gate→ turbine.opening
#
# Output files:
#   images/connection_test_v2.png     ← time-domain response
#   images/frequency_response_v2.png  ← Bode: S(s) and T(s)
#   images/nyquist_v2.png             ← Nyquist: L(s) with stability margins
# =============================================================================

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks           # get_sensitivity, get_looptransfer
using OrdinaryDiffEq
using ControlSystemsBase                               # ss(), margin(), bodeplot()
using Plots
using Printf

# ── Operating Point B (Trollheim benchmark) ───────────────────────────────────
const H_n       = 371.0
const Q_n       = 37.0
const P_r       = 130.0e6
const n_rpm     = 375.0
const omega_ref = 2π * n_rpm / 60    # ≈ 39.27 rad/s
const ρ         = 1000.0
const g_acc     = 9.81

# ── Penstock geometry (lumped two-section equivalent) ─────────────────────────
const L_p = 363.0 + 145.0   # 508 m total
const D_p = 4.1              # area-weighted equivalent diameter [m]

# ── Turbine sizing ────────────────────────────────────────────────────────────
const D_runner  = 2.5
const τ_o_ss    = 0.90
const K_q_v     = Q_n / (τ_o_ss * D_runner^2 * sqrt(H_n))   # ≈ 0.3415
const η_max_v   = 0.97
const c_eta_v   = 0.25
const Q_rated_v = Q_n

# ── Rotor (M_a = 8 s) ────────────────────────────────────────────────────────
const J_rotor = 8.0 * P_r / omega_ref^2    # ≈ 6.74×10⁵ kg·m²

# ── GGOV1 droop (σ = 0.10) ───────────────────────────────────────────────────
const sigma       = 0.10
const T_ps        = 1.75
const T_gs        = 0.20
const K_i_gov     = 1.0 / T_gs               # gate integrator gain = 5.0
const R_droop_phy = sigma * omega_ref / P_r  # ≈ 3.02e-8 [rad/s/W]

println("=== Trollheim HPP — Component Model v2 (Sensors + Analysis Points) ===")
@printf "  ω_ref   = %.4f rad/s  (%g RPM)\n"   omega_ref n_rpm
@printf "  K_q     = %.4f\n"                    K_q_v
@printf "  J_rotor = %.3e kg·m²\n"              J_rotor
@printf "  sigma   = %.2f   T_gov = %.2f s   K_i = %.2f\n" sigma T_ps K_i_gov
println()

# =============================================================================
# ── Composite model (MTKSL-style: sensors + analysis points) ─────────────────
# =============================================================================
@mtkmodel TrollheimDroopV2 begin
    @components begin
        # ── Hydraulic circuit ─────────────────────────────────────────────────
        reservoir  = Reservoir(H = H_n, rho = ρ, g = g_acc)
        penstock   = Penstock(L = L_p, D_pipe = D_p, rho = ρ, e = 1.5e-5)
        turbine    = FrancisTurbineAffinity(
                         D = D_runner, K_q = K_q_v,
                         eta_max = η_max_v, c_eta = c_eta_v,
                         Q_rated = Q_rated_v, rho = ρ, g = g_acc)
        tailwater  = Reservoir(H = 0.0, rho = ρ, g = g_acc)

        # ── Mechanical train ──────────────────────────────────────────────────
        rotor      = RotorInertia(J = J_rotor, b_v = 50.0, omega_0 = omega_ref)
        generator  = SimpleGenerator(
                         P_rated = P_r, omega_rated = omega_ref,
                         omega_s = omega_ref, D_d = 1.5, eta_gen = 0.985)

        # ── Speed sensor (RotationalPort → SignalOutPort) ─────────────────────
        # Attached in parallel on the flange_generator node.
        # Draws zero torque (ideal probe) — does not disturb mechanical balance.
        speed_sensor = RotationalSpeedSensor()

        # ── GGOV1 droop governor ──────────────────────────────────────────────
        governor   = GGOV1Governor(
                         R_droop = R_droop_phy, T_p = 0.05,
                         T_gov = T_ps, omega_ref = omega_ref,
                         P_set = P_r, tau_min = 0.0, tau_max = 1.0,
                         K_i = K_i_gov, tau_0 = τ_o_ss)
    end

    @equations begin
        # ── Hydraulic (acausal, Kirchhoff-type) ───────────────────────────────
        connect(reservoir.port,   penstock.port_a)
        connect(penstock.port_b,  turbine.port_in)
        connect(turbine.port_out, tailwater.port)

        # ── Mechanical (acausal, Kirchhoff-type) ──────────────────────────────
        connect(turbine.shaft,         rotor.flange_turbine)

        # 3-way connect: speed_sensor taps in parallel at flange_generator node.
        # MTK Kirchhoff: sum(tau) = 0  →  speed_sensor contributes 0 (ideal).
        # Both generator and speed_sensor receive same omega as rotor.
        connect(rotor.flange_generator, speed_sensor.flange, generator.flange)

        # ── Cross-domain with analysis points (causal signal) ─────────────────
        #
        # :y_omega — plant OUTPUT — measured speed signal leaving the physical
        #            domain and entering the governor's speed input.
        #            This is the "y" in the classic control feedback loop.
        connect(speed_sensor.w, :y_omega, governor.speed_in)

        # :u_gate  — plant INPUT — control action leaving the governor and
        #            entering the turbine's hydraulic guide-vane.
        #            This is the "u" in the classic control feedback loop.
        connect(governor.gate_out, :u_gate, turbine.opening)

        # Power feedback (inner measurement — not the primary analysis loop,
        # kept as a simple algebraic connection, no analysis point needed).
        governor.power_in.u ~ turbine.P_mech
    end
end

# ── Build named (uncompiled) model — MUST come before mtkcompile ──────────────
# The uncompiled `model` is what Blocks.get_sensitivity needs to locate :y_omega
@named model = TrollheimDroopV2()

# ── Compile ───────────────────────────────────────────────────────────────────
sys = mtkcompile(model)

println("  States after compilation: ", length(unknowns(sys)))
println()

# ── Steady-state initial conditions ──────────────────────────────────────────
Q_ic    = τ_o_ss * K_q_v * D_runner^2 * sqrt(H_n)
P_ic    = ρ * g_acc * Q_ic * H_n * η_max_v
dm_ic   = ρ * Q_ic
η_ic    = η_max_v * (1.0 - c_eta_v * (Q_ic / Q_rated_v - 1.0)^2)
A_pipe  = π * D_p^2 / 4
μ_water = 1e-3
Re_ic   = dm_ic * D_p / (μ_water * A_pipe)
f_D_ic  = 1.0 / (2.0 * log10(1.5e-5 / (3.7 * D_p) + 5.7 / Re_ic^0.9))^2
dp_f_ic = f_D_ic * (L_p / D_p) * (dm_ic / (ρ * A_pipe))^2 * ρ / 2

u0 = Dict(
    sys.rotor.omega      => omega_ref,
    sys.rotor.theta      => 0.0,
    sys.penstock.dm      => dm_ic,
    sys.penstock.p_avg   => ρ * g_acc * H_n,
    # MTSL block states (pow_lag / gov_lag are FirstOrder, gate_int is Integrator).
    # The alias variables P_meas / x_gov / x_int are now algebraic — use .x directly.
    sys.governor.pow_lag.x   => P_ic,    # pow_lag output ≡ P_meas
    sys.governor.gov_lag.x   => 0.0,     # gov_lag output ≡ x_gov
    sys.governor.gate_int.x  => τ_o_ss,  # integrator state ≡ x_int
)

guesses = Dict(
    sys.turbine.H         => H_n,
    sys.turbine.Q         => Q_ic,
    sys.turbine.eta       => η_ic,
    sys.turbine.P_mech    => P_ic,
    sys.turbine.dm        => dm_ic,
    sys.turbine.tau_shaft => P_ic / omega_ref,
    sys.penstock.Re       => Re_ic,
    sys.penstock.f_D      => f_D_ic,
    sys.penstock.dp_f     => dp_f_ic,
)

tspan = (0.0, 60.0)
prob  = ODEProblem(sys, u0, tspan; guesses)

# ── 10% load step at t = 5 s ─────────────────────────────────────────────────
cb = DiscreteCallback(
    (u, t, i) -> t == 5.0,
    i -> (i.ps[sys.governor.P_set] = 1.1 * P_r; nothing),
    save_positions = (true, true))

# ── Solve ─────────────────────────────────────────────────────────────────────
println("Solving time-domain simulation (0 → 60 s)...")
@time sol = solve(prob, Rodas5P(); tstops = [5.0], callback = cb,
                  abstol = 1e-7, reltol = 1e-7)
println("  Retcode: ", sol.retcode)

ω_ss   = sol(60.0; idxs = sys.rotor.omega)
Δω_ss  = ω_ss - omega_ref
Δω_pred = -0.1 * sigma * omega_ref
@printf "  ω_ss (t=60)      = %.5f rad/s   (ref = %.5f)\n"    ω_ss omega_ref
@printf "  Δω_ss measured   = %+.5f rad/s\n"                  Δω_ss
@printf "  Δω_ss predicted  = %+.5f rad/s  (σ=%.2f, 10%% step)\n" Δω_pred sigma
println()

# ── Time-domain plots ─────────────────────────────────────────────────────────
default(linewidth = 1.8, tickfontsize = 9, guidefontsize = 10,
        legendfontsize = 8, left_margin = 8Plots.mm, bottom_margin = 3Plots.mm)

tv   = range(0.0, 60.0, 600)
ω_v  = [sol(τ; idxs = sys.rotor.omega)          for τ in tv]
τo_v = [sol(τ; idxs = sys.governor.tau_o)       for τ in tv]
Pm_v = [sol(τ; idxs = sys.turbine.P_mech) / 1e6 for τ in tv]
Δω_v = [(ω - omega_ref) / omega_ref * 100        for ω in ω_v]

p_ω  = plot(tv, ω_v;  label = "ω [rad/s]",     color = :steelblue,
            ylabel = "Speed  [rad/s]",
            title  = "Trollheim HPP v2 — GGOV1 Droop (σ=0.10) | 10% step at t=5 s")
hline!([omega_ref];           ls = :dash, lc = :black,  lw = 1, label = "ω_ref")
hline!([omega_ref + Δω_pred]; ls = :dot,  lc = :orange, lw = 1, label = "Δω_ss predicted")
vline!([5.0];                 ls = :dot,  lc = :grey,   lw = 1, label = "step")

p_g  = plot(tv, τo_v; label = "τ_o [-]", color = :steelblue,
            ylabel = "Gate  [pu]", legend = :bottomright)
vline!([5.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p_P  = plot(tv, Pm_v; label = "P_mech [MW]", color = :steelblue,
            ylabel = "Mech. power  [MW]", legend = :topleft)
hline!([P_r / 1e6];       ls = :dash, lc = :black,  lw = 1, label = "P_rated")
hline!([1.1 * P_r / 1e6]; ls = :dot,  lc = :orange, lw = 1, label = "1.1 P_r")
vline!([5.0]; ls = :dot, lc = :grey, lw = 1, label = false)

p_dω = plot(tv, Δω_v; label = "Δω/ω_ref [%]", color = :steelblue,
            ylabel = "Freq. deviation  [%]", legend = :bottomleft)
hline!([Δω_pred / omega_ref * 100]; ls = :dot,  lc = :orange, lw = 1, label = "predicted SS")
hline!([0.0];                       ls = :dash, lc = :black,  lw = 1, label = "0 %")
vline!([5.0]; ls = :dot, lc = :grey, lw = 1, label = false)

mkpath(joinpath(@__DIR__, "images"))
savefig(plot(p_ω, p_g, p_P, p_dω; layout = (4,1), size = (950,1100),
             link = :x, xlabel = "Time  [s]"),
        joinpath(@__DIR__, "images", "connection_test_v2.png"))
println("Time-domain plot saved  →  images/connection_test_v2.png")
println()

# =============================================================================
# ── FREQUENCY-DOMAIN ANALYSIS via ControlSystemsBase ─────────────────────────
#
# Linearise around the steady-state operating point using the analysis points:
#   :y_omega — plant output  (speed measurement leaving physical domain)
#   :u_gate  — plant input   (gate command entering turbine)
#
# The linearisation point is the pre-step steady state (gate=0.90, rated speed).
# We build `op` by merging the known IC values into a zero-padded dict of all
# compiled unknowns.  MTK evaluates the Jacobian of the simplified equations at
# this point to obtain the state-space representation A, B, C, D.
# =============================================================================
println("Linearising at steady-state operating point...")

# Operating point: compiled-state dict, padded with zeros for unspecified vars
op_known = Dict(
    sys.rotor.omega     => omega_ref,
    sys.rotor.theta     => 0.0,
    sys.penstock.dm     => dm_ic,
    sys.penstock.p_avg  => ρ * g_acc * H_n,
    sys.governor.pow_lag.x  => P_ic,
    sys.governor.gov_lag.x  => 0.0,
    sys.governor.gate_int.x => τ_o_ss,
)
op = merge(Dict(unknowns(sys) .=> 0.0), op_known)

# ── Sensitivity  S(s) = 1 / (1 + P·C) ───────────────────────────────────────
# This captures how well disturbances at the plant output are rejected.
# A peak |S(jω)| > 6 dB signals poor disturbance rejection.
matrices_S, _ = Blocks.get_sensitivity(model, :y_omega; op)
S_sys = ss(matrices_S...) |> minreal

# ── Complementary sensitivity  T(s) = P·C / (1 + P·C) ───────────────────────
# Reference tracking quality.  T(0) → 1 for integral governor (isochronous).
matrices_T, _ = Blocks.get_comp_sensitivity(model, :y_omega; op)
T_sys = ss(matrices_T...)

# ── Loop transfer function  L(s) = P·C ───────────────────────────────────────
# Note: ControlSystemsBase convention for negative-feedback Nyquist uses -L.
matrices_L, _ = Blocks.get_looptransfer(model, :y_omega; op)
L_sys = -ss(matrices_L...)   # negate built-in negative feedback sign

# ── Stability margins ─────────────────────────────────────────────────────────
gm_lin, pm_lin, ωgc, ωpc = margin(L_sys)
@printf "  Gain  margin  GM = %.2f dB   at ω = %.4f rad/s\n"  20log10(gm_lin) ωpc
@printf "  Phase margin  PM = %.2f °    at ω = %.4f rad/s\n"  pm_lin          ωgc
@printf "  (targets: GM > 6 dB,  PM > 45°)\n"
println()

Ms, ωMs = hinfnorm(S_sys)   # peak sensitivity (robustness measure)
@printf "  Peak sensitivity  |S|_∞ = %.3f  (target < 6 dB = 2.0)\n" Ms
println()

# ── Bode plot: S(s) and T(s) ─────────────────────────────────────────────────
p_bode = bodeplot([S_sys, T_sys];
                  label     = ["S(s)  sensitivity" "T(s)  comp. sensitivity"],
                  plotphase = false,
                  title     = "Trollheim GGOV1  σ=$(sigma) — Sensitivity functions",
                  yscale    = :identity)   # magnitude in dB by default

savefig(p_bode, joinpath(@__DIR__, "images", "frequency_response_v2.png"))
println("Bode plot saved  →  images/frequency_response_v2.png")

# ── Nyquist plot of L(s) ──────────────────────────────────────────────────────
p_nyq = nyquistplot(L_sys;
                    label      = "L(s)",
                    title      = "Trollheim GGOV1  σ=$(sigma) — Nyquist diagram",
                    Ms_circles = Ms)   # draw robustness circle at peak sensitivity

savefig(p_nyq, joinpath(@__DIR__, "images", "nyquist_v2.png"))
println("Nyquist plot saved  →  images/nyquist_v2.png")
println()
println("Done.")
