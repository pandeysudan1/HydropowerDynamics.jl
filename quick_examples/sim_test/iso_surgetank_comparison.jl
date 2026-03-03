# =============================================================================
# iso_surgetank_comparison.jl
# Trollheim HPP — Effect of Surge Tank on Isochronous Governor Response
#
# Two scenarios, both with an ISO PID governor.
# Gate position ramp: tau_o  0.50 pu  →  0.95 pu  (triggered by load step at t = 50 s)
#
#   Scenario 1  —  ISO without surge tank
#                  Constant net head H = H_n = 371 m  (penstock only)
#
#   Scenario 2  —  ISO with surge tank
#                  Headrace tunnel + surge tank + penstock:
#                    dQ_t/dt : headrace momentum (water inertia)
#                    dZ/dt   : surge tank free-surface (SurgeTank  ST.1–ST.2)
#                  Net head at turbine oscillates at the surge natural period
#                    T_s = 2π√(L_t · A_surge / (g · A_tunnel)) ≈ 73 s
#
# Waterway geometry from:
#   benchmarks/trollheim_hydropower.md  (Energies 2019 table)
#
# Output:
#   images/iso_surgetank_gate_ramp.png
# =============================================================================

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots
using Printf

# ── Operating Point B ────────────────────────────────────────────────────────
const H_n     = 371.0              # nominal net head        [m]
const Q_n     = 37.0               # nominal discharge       [m³/s]
const n_rpm   = 375.0              # nominal speed           [RPM]
const ω_ref   = 2π * n_rpm / 60   # ≈ 39.27 rad/s
const ρ       = 1000.0             # water density           [kg/m³]
const g_acc   = 9.81               # gravity                 [m/s²]

# ── Turbine parameters ───────────────────────────────────────────────────────
const D_r       = 2.5                          # runner diameter       [m]
const η_max_v   = 0.97                         # peak efficiency       [-]
const c_eta_v   = 0.25                         # hill curve width      [-]
const Q_rated_v = Q_n                          # rated flow at BEP     [m³/s]

# K_q sized so that tau_o=0.90 gives Q_n at H_n  (design point, unchanged)
const τ_o_design = 0.90
const K_q_v      = Q_n / (τ_o_design * D_r^2 * sqrt(H_n))  # ≈ 0.3415

# ── Initial operating point: gate = 0.50 ─────────────────────────────────────
const τ_o_ss    = 0.50
const Q_ss      = τ_o_ss * K_q_v * D_r^2 * sqrt(H_n)        # ≈ 20.5 m³/s
const η_ss      = η_max_v * (1.0 - c_eta_v * (Q_ss/Q_rated_v - 1.0)^2)
const P_mech_ss = ρ * g_acc * Q_ss * H_n * η_ss              # ≈ 69 MW

# ── Target operating point: gate = 0.95 ──────────────────────────────────────
const τ_o_step  = 0.95
const Q_step    = τ_o_step * K_q_v * D_r^2 * sqrt(H_n)       # ≈ 39.1 m³/s
const η_step    = η_max_v * (1.0 - c_eta_v * (Q_step/Q_rated_v - 1.0)^2)
const P_mech_step = ρ * g_acc * Q_step * H_n * η_step         # ≈ 138 MW

# ── Rotor & generator ────────────────────────────────────────────────────────
# Use rated J (based on full-load power) so inertia is physically correct
const P_rated_design = ρ * g_acc * Q_n * H_n * η_max_v   # ≈ 130.6 MW
const J_v     = 8.0 * P_rated_design / ω_ref^2           # M_a = 8 s → ≈ 6.74×10⁵ kg·m²

# ── Governor (ISO PID, Trollheim tuning) ────────────────────────────────────
const Kp          = 2.0
const Ki          = 2.0
const e_int_ss    = τ_o_ss / Ki        # integrator IC so u_pid = τ_o_ss at t=0 (= 0.25)
const R_open_v    = 0.05               # max gate open  rate  [pu/s]
const R_close_v   = 0.20               # max gate close rate  [pu/s]

# ── Waterway geometry  (Trollheim benchmark table) ───────────────────────────
# Headrace tunnel: three sections → total length 4 496.5 m, diameter 6.3 m
const L_tunnel   = 81.5 + 395.0 + 4020.0          # 4 496.5 m
const D_tunnel   = 6.3                             # m
const A_tunnel   = π/4 * D_tunnel^2               # 31.17 m²

# Surge tank: diameter 3.4 m, vertical extent 75.5 m
const D_st       = 3.4                             # surge tank diameter [m]
const A_surge    = π/4 * D_st^2                   # 9.08 m²

# Penstock: two sections (from waterway table)
#   #1: L = 363 m, D = 4.7 m   #2: L = 145 m, D = 3.3 m
const A_p1       = π/4 * 4.7^2                    # 17.35 m²
const A_p2       = π/4 * 3.3^2                    # 8.55 m²
const R_p_v      = 0.015 * 363.0 / (4.7 * 2 * g_acc * A_p1^2) +
                   0.015 * 145.0 / (3.3 * 2 * g_acc * A_p2^2)  # ≈ 6.55×10⁻⁴ s²/m⁵
const R_t_v      = 0.012 * L_tunnel / (D_tunnel * 2 * g_acc * A_tunnel^2)  # ≈ 4.49×10⁻⁴ s²/m⁵

# ── Surge tank steady-state level (reference: tailwater = 0 m) ──────────────
# At SS (gate=0.50): H_eff = Z_0 - R_p·Q_ss² = H_n  →  Z_0 = H_n + R_p·Q_ss²
const Z_0        = H_n + R_p_v * Q_ss^2          # ≈ 371.2 m  (initial SS, gate=0.50)
const H_res      = Z_0 + R_t_v * Q_ss^2          # ≈ 371.5 m (reservoir equiv. head)

# ── Surge natural period (analytical) ────────────────────────────────────────
const T_surge    = 2π * sqrt(L_tunnel * A_surge / (g_acc * A_tunnel))  # ≈ 73 s

println("=== Trollheim — ISO Governor + Surge Tank Comparison ===")
@printf "  Gate ramp        : tau_o  %.2f  ->  %.2f pu  (step at t=50 s)\n"  τ_o_ss  τ_o_step
@printf "  ω_ref            = %.4f rad/s  (%g RPM)\n"  ω_ref   n_rpm
@printf "  K_q              = %.4f\n"                      K_q_v
@printf "  Q_ss  (gate=%.2f) = %.2f m³/s  →  P_mech = %.1f MW\n"  τ_o_ss  Q_ss  P_mech_ss/1e6
@printf "  Q_step(gate=%.2f) = %.2f m³/s  →  P_mech = %.1f MW\n"  τ_o_step Q_step P_mech_step/1e6
@printf "  J_rotor          = %.3e  kg·m²\n"               J_v
println()
@printf "  L_tunnel         = %.1f m     A_tunnel = %.2f m²\n" L_tunnel A_tunnel
@printf "  A_surge_tank     = %.2f m²   (D = %.1f m)\n"       A_surge  D_st
@printf "  R_t              = %.3e s²/m⁵\n"                   R_t_v
@printf "  R_p              = %.3e s²/m⁵\n"                   R_p_v
@printf "  Z_0 (SS level)   = %.2f m    H_res = %.2f m\n"     Z_0  H_res
@printf "  Surge period T_s = %.1f s\n"                       T_surge
println()

# =============================================================================
# ── Scenario 1 : ISO — no surge tank  (constant H = H_n) ─────────────────────
# =============================================================================
@mtkmodel TrollheimISONO begin
    @parameters begin
        H        = H_n
        rho      = ρ;       g        = g_acc
        D_runner = D_r;     K_q      = K_q_v
        eta_max  = η_max_v; c_eta    = c_eta_v;  Q_rated = Q_rated_v
        J        = J_v
        P_rated  = P_mech_ss        # stepped +10 % by callback
        K_p      = Kp;      K_i     = Ki
        omega_ref= ω_ref
        R_open   = R_open_v;  R_close = R_close_v
    end
    @variables begin
        omega(t)     = ω_ref
        Q(t)
        eta(t)
        P_mech(t)
        tau_turb(t)
        tau_gen(t)
        e(t)
        e_int(t)     = e_int_ss
        tau_o(t)     = τ_o_ss
        u_pid(t)
    end
    @equations begin
        # FA.1–FA.3 (constant head)
        Q       ~  tau_o * K_q * D_runner^2 * sqrt(abs(H) + 1e-6)
        eta     ~  eta_max * (1.0 - c_eta * (Q / Q_rated - 1.0)^2)
        P_mech  ~  rho * g * Q * H * eta
        tau_turb~  P_mech / (abs(omega) + 1e-6)
        # Generator (GA.1)
        tau_gen ~  P_rated / omega_ref
        # Swing equation (M.1)
        D(omega)~  (tau_turb - tau_gen) / J
        # ISO PID governor
        e       ~  omega_ref - omega
        D(e_int)~  e
        u_pid   ~  K_p * e + K_i * e_int
        D(tau_o)~  max(-R_close, min(R_open,
                       (clamp(u_pid, 0.0, 1.0) - tau_o) / 0.01))
    end
end

# =============================================================================
# ── Scenario 2 : ISO — WITH surge tank  ──────────────────────────────────────
#
# Added states:   Q_t [m³/s]  headrace tunnel volumetric flow
#                 Z   [m]     surge tank free-surface level
#
# Surge-tank equations  (SurgeTank  ST.1–ST.2, hydraulic.jl §2.4):
#   dZ/dt     = (Q_t - Q) / A_surge                          (ST.1)
#   H_eff     = Z  (ST.2 with p_atm cancelling at turbine)
#
# Penstock algebraic solution (quasi-static, valid for T_surge >> T_w):
#   FA.1 + H_eff = Z - R_p·Q²  →  H_eff = Z/(1 + R_p·φ²)
#   where  φ = tau_o·K_q·D²    (analytical elimination of algebraic loop)
#
# Headrace momentum  (Penstock P.1, §2.3, low-friction tunnel):
#   dQ_t/dt   = (g·A_tunnel/L_tunnel)·(H_res - Z - R_t·Q_t·|Q_t|)
# =============================================================================
@mtkmodel TrollheimISOSurge begin
    @parameters begin
        rho      = ρ;       g        = g_acc
        D_runner = D_r;     K_q      = K_q_v
        eta_max  = η_max_v; c_eta    = c_eta_v;  Q_rated = Q_rated_v
        J        = J_v
        P_rated  = P_mech_ss
        K_p      = Kp;      K_i     = Ki
        omega_ref= ω_ref
        R_open   = R_open_v;  R_close = R_close_v
        # Waterway surge parameters
        H_reservoir = H_res
        R_t      = R_t_v         # headrace friction   [s²/m⁵]
        R_p      = R_p_v         # penstock friction   [s²/m⁵]
        A_t      = A_tunnel      # tunnel cross-section [m²]
        L_t      = L_tunnel      # tunnel length        [m]
        A_s      = A_surge       # surge tank x-sect    [m²]
    end
    @variables begin
        omega(t)     = ω_ref
        Q(t)
        H_eff(t)
        eta(t)
        P_mech(t)
        tau_turb(t)
        tau_gen(t)
        e(t)
        e_int(t)     = e_int_ss
        tau_o(t)     = τ_o_ss
        u_pid(t)
        Q_t(t)       = Q_ss         # headrace tunnel flow IC = SS flow at gate=0.50 [m³/s]
        Z(t)         = Z_0          # surge tank water level [m]
    end
    @equations begin
        # ── Penstock quasi-static algebraic solution (ST.2 + FA.1 combined)
        # phi = tau_o * K_q * D²  →  H_eff = Z / (1 + R_p * phi^2)
        H_eff  ~  Z / (1.0 + R_p * (tau_o * K_q * D_runner^2)^2 + 1e-10)

        # ── FA.1–FA.3 (variable head H_eff from surge level Z)
        Q       ~  tau_o * K_q * D_runner^2 * sqrt(abs(H_eff) + 1e-6)
        eta     ~  eta_max * (1.0 - c_eta * (Q / Q_rated - 1.0)^2)
        P_mech  ~  rho * g * Q * H_eff * eta
        tau_turb~  P_mech / (abs(omega) + 1e-6)

        # ── Generator
        tau_gen ~  P_rated / omega_ref

        # ── Swing equation (M.1)
        D(omega)~  (tau_turb - tau_gen) / J

        # ── ISO PID governor
        e       ~  omega_ref - omega
        D(e_int)~  e
        u_pid   ~  K_p * e + K_i * e_int
        D(tau_o)~  max(-R_close, min(R_open,
                       (clamp(u_pid, 0.0, 1.0) - tau_o) / 0.01))

        # ── Surge tank free-surface  (ST.1, hydraulic.jl §2.4)
        D(Z)    ~  (Q_t - Q) / A_s

        # ── Headrace tunnel momentum  (Penstock P.1, §2.3)
        D(Q_t)  ~  (g * A_t / L_t) *
                   (H_reservoir - Z - R_t * Q_t * sqrt(Q_t^2 + 1e-6))
    end
end

# ── Compile both models ───────────────────────────────────────────────────────
@mtkcompile sys1 = TrollheimISONO()
@mtkcompile sys2 = TrollheimISOSurge()

# ── Initial conditions and parameters ────────────────────────────────────────
u0_s1 = [sys1.omega => ω_ref,  sys1.e_int => e_int_ss,  sys1.tau_o => τ_o_ss]
u0_s2 = [sys2.omega => ω_ref,  sys2.e_int => e_int_ss,  sys2.tau_o => τ_o_ss,
          sys2.Q_t  => Q_ss,    sys2.Z    => Z_0]    # Q_t IC matches gate=0.50 SS flow

p0_s1 = [sys1.P_rated => P_mech_ss]
p0_s2 = [sys2.P_rated => P_mech_ss]

tspan = (0.0, 400.0)   # 400 s — covers ≈ 5 surge oscillation periods

prob1 = ODEProblem(sys1, merge(Dict(u0_s1), Dict(p0_s1)), tspan)
prob2 = ODEProblem(sys2, merge(Dict(u0_s2), Dict(p0_s2)), tspan)

# ── Load step at t = 50 s: P_rated jumps from P_mech_ss (gate=0.50) ──────────
#    to P_mech_step (gate=0.95), forcing governor to open gate fully
cb1 = DiscreteCallback((u,t,i) -> t == 50.0,
        i -> (i.ps[sys1.P_rated] = P_mech_step; nothing),
        save_positions = (true,true))

cb2 = DiscreteCallback((u,t,i) -> t == 50.0,
        i -> (i.ps[sys2.P_rated] = P_mech_step; nothing),
        save_positions = (true,true))

# ── Solve ─────────────────────────────────────────────────────────────────────
println("Solving Scenario 1 (no surge tank)...")
@time sol1 = solve(prob1, Rodas5P(), tstops=[50.0], callback=cb1,
                   abstol=1e-7, reltol=1e-7)
println("  Retcode: ", sol1.retcode)
@printf "  tau_o_ss (t=400) = %.4f pu  (target=%.2f)\n"  sol1[sys1.tau_o][end]  τ_o_step
@printf "  ω_ss     (t=400) = %.5f  (ref=%.5f,  Δω=%.5e rad/s)\n" sol1[sys1.omega][end] ω_ref (sol1[sys1.omega][end]-ω_ref)

println()
println("Solving Scenario 2 (with surge tank)...")
@time sol2 = solve(prob2, Rodas5P(), tstops=[50.0], callback=cb2,
                   abstol=1e-7, reltol=1e-7)
println("  Retcode: ", sol2.retcode)
@printf "  tau_o_ss (t=400) = %.4f pu  (target=%.2f)\n"  sol2[sys2.tau_o][end]  τ_o_step
@printf "  ω_ss     (t=400) = %.5f  (ref=%.5f,  Δω=%.5e rad/s)\n" sol2[sys2.omega][end] ω_ref (sol2[sys2.omega][end]-ω_ref)
@printf "  Z range: %.2f - %.2f m  (surge oscillation amplitude = %.3f m)\n" minimum(sol2[sys2.Z]) maximum(sol2[sys2.Z]) (maximum(sol2[sys2.Z])-minimum(sol2[sys2.Z]))/2
println()

# ── Extract time series ───────────────────────────────────────────────────────
tv1  = sol1.t;   tv2  = sol2.t
ω1   = sol1[sys1.omega];       ω2   = sol2[sys2.omega]
τo1  = sol1[sys1.tau_o];       τo2  = sol2[sys2.tau_o]
Q1   = sol1[sys1.Q];           Q2   = sol2[sys2.Q]
Pm1  = sol1[sys1.P_mech]./1e6; Pm2  = sol2[sys2.P_mech]./1e6
Z2   = sol2[sys2.Z]
H1   = fill(H_n, length(tv1)); H2   = sol2[sys2.H_eff]

# ── Plot ──────────────────────────────────────────────────────────────────────
default(linewidth = 1.8, tickfontsize = 9, guidefontsize = 10,
        legendfontsize = 8, left_margin = 8Plots.mm, bottom_margin = 3Plots.mm)

c1 = :steelblue;   c2 = :firebrick

p_omega = plot(tv1, ω1; label="S1: no surge tank", color=c1,
               ylabel="Speed  [rad/s]",
               title="ISO Governor — Surge Tank Effect  |  Trollheim HPP  (gate $(τ_o_ss) -> $(τ_o_step) pu  at t=50 s)")
plot!(p_omega, tv2, ω2; label="S2: with surge tank", color=c2)
hline!([ω_ref]; ls=:dash, lc=:black, lw=1, label="omega_ref")
vline!([50.0];  ls=:dot,  lc=:grey,  lw=1, label="step")

p_gate = plot(tv1, τo1; label="S1: no surge tank", color=c1,
              ylabel="Gate  [pu]", legend=:bottomright)
plot!(p_gate, tv2, τo2; label="S2: with surge tank", color=c2)
vline!([50.0]; ls=:dot, lc=:grey, lw=1, label=false)

p_head = plot(tv1, H1; label="S1: H = const = $(H_n) m", color=c1, ls=:dash,
              ylabel="Net head  [m]", legend=:topleft)
plot!(p_head, tv2, H2; label="S2: H_eff (surge oscillation)", color=c2)
vline!([50.0]; ls=:dot, lc=:grey, lw=1, label=false)

p_q = plot(tv1, Q1; label="S1: no surge tank", color=c1,
           ylabel="Discharge  [m3/s]", legend=:topleft)
plot!(p_q, tv2, Q2; label="S2: with surge tank", color=c2)
vline!([50.0]; ls=:dot, lc=:grey, lw=1, label=false)

p_z = plot(tv2, Z2; label="Surge level Z(t)", color=:darkgreen,
           ylabel="Surge level  [m]", legend=:topleft)
hline!([Z_0]; ls=:dash, lc=:black, lw=1, label="Z_0 = $(round(Z_0,digits=1)) m  (SS)")
vline!([50.0]; ls=:dot, lc=:grey, lw=1, label="step")

ptop = plot(p_omega, p_gate, p_head, p_q, p_z;
            layout       = (5, 1),
            size         = (950, 1300),
            link         = :x,
            xlabel       = "Time  [s]")

mkpath(joinpath(@__DIR__, "images"))
out_png = joinpath(@__DIR__, "images", "iso_surgetank_gate_ramp.png")
savefig(ptop, out_png)
println("Plot saved  ->  sim_test/images/iso_surgetank_gate_ramp.png")
