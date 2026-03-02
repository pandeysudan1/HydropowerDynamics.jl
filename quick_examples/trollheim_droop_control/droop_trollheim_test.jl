# =============================================================================
# droop_trollheim_test.jl
# Trollheim HPP — Droop (GGOV1-style) Governor Step Response
#
# Operating Point B:  H = 371 m,  Q = 37 m³/s,  P ≈ 130 MW,  n = 375 RPM
# Scenario:           10 % load increase applied at t = 50 s
# Simulation:         0 s → 100 s   (steady state before t = 50 s)
#
# Governor type:      GGOV1-style with droop σ = 0.1 (Trollheim benchmark)
# Expected outcome:   speed settles at ω_ref − Δω_ss  (≈ −1 % deviation)
#
# Benchmark data (Trollheim, governor table):
#   σ = 0.10      T_w  = 0.44 s    T_gs = 0.20 s    T_ps = 1.75 s
#   R_open = 0.05 pu/s             R_close = 0.20 pu/s
# =============================================================================

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots
using Printf

mkpath(joinpath(@__DIR__, "quick_examples", "images"))

# ── Plant constants (same as ISO test) ───────────────────────────────────────
const H_n_d    = 371.0
const Q_n_d    = 37.0
const n_rpm_d  = 375.0
const ω_ref_d  = 2π * n_rpm_d / 60   # ≈ 39.27 rad/s
const ρ_d      = 1000.0
const g_d      = 9.81

# Turbine
const D_r_d     = 2.5
const η_max_d   = 0.97
const c_eta_d   = 0.25
const Q_rated_d = Q_n_d
const τ_o_ss_d  = 0.90
const K_q_d     = Q_n_d / (τ_o_ss_d * D_r_d^2 * sqrt(H_n_d))
const P_mech_ss_d = ρ_d * g_d * Q_n_d * H_n_d * η_max_d   # ≈ 130.6 MW
const J_d       = 8.0 * P_mech_ss_d / ω_ref_d^2

# Droop governor parameters (Trollheim benchmark table)
const σ_d    = 0.10    # static droop (regulation)
const T_p_d  = 0.05    # power measurement lag [s]
const T_gov_d = 1.75   # pilot servo / governor lag = T_ps [s]
const K_i_d  = 1.0 / 0.20   # gate integrator gain = 1/T_gs = 5.0 [1/s]

# Droop coefficient in physical units:  ΔP [W] = P_r/(σ·ω_ref) · Δω [rad/s]
# → coefficient K_droop = P_mech_ss / (σ * ω_ref)  [W·s/rad]
const K_droop_d = P_mech_ss_d / (σ_d * ω_ref_d)

# Predicted steady-state speed deviation after 10% load step
# At new SS:  K_droop * Δω_ss = ΔP_load = 0.1 * P_mech_ss
# → Δω_ss = 0.1 * P_mech_ss / K_droop = 0.1 * σ * ω_ref = −0.01·ω_ref
const Δω_predicted = -(0.10 * σ_d * ω_ref_d)  # ≈ −0.393 rad/s  (−1 %)

println("=== Trollheim Droop Governor ===")
@printf "  ω_ref           = %.4f rad/s  (%g RPM)\n" ω_ref_d n_rpm_d
@printf "  K_q             = %.4f\n" K_q_d
@printf "  P_mech_ss       = %.2f MW\n" P_mech_ss_d/1e6
@printf "  J_rotor         = %.0f kg·m²\n" J_d
@printf "  σ (droop)       = %.2f\n" σ_d
@printf "  T_gov (T_ps)    = %.2f s\n" T_gov_d
@printf "  Predicted Δω_ss = %.4f rad/s  (%.2f %%)\n" Δω_predicted Δω_predicted/ω_ref_d*100
println()

# ── Flat MTK model ────────────────────────────────────────────────────────────
@mtkmodel TrollheimDroop begin
    @parameters begin
        # Hydraulic / turbine
        H        = H_n_d
        rho      = ρ_d
        g        = g_d
        D_runner = D_r_d
        K_q      = K_q_d
        eta_max  = η_max_d
        c_eta    = c_eta_d
        Q_rated  = Q_rated_d
        # Rotor
        J        = J_d
        D_d_gen  = 2.0           # generator synchronising torque coeff
        # Generator load (stepped at t=50 via callback)
        P_load   = P_mech_ss_d   # generator electrical load [W] – STEPPED
        # Droop governor – setpoint is FIXED (NOT stepped)
        P_set    = P_mech_ss_d   # governor power setpoint [W] – constant
        K_droop  = K_droop_d     # physical droop coeff [W·s/rad]
        T_p      = T_p_d
        T_gov    = T_gov_d
        K_i_gov  = K_i_d
        omega_ref = ω_ref_d
    end
    @variables begin
        omega(t)   = ω_ref_d       # rotor speed [rad/s]
        Q(t)                       # flow [m³/s]
        eta(t)                     # efficiency [-]
        P_mech(t)                  # mechanical power [W]
        tau_turbine(t)             # turbine torque [N·m]
        tau_gen(t)                 # generator reaction torque [N·m]
        # GGOV1 governor states
        P_ref(t)                   # load reference with droop [W]
        P_meas(t) = P_mech_ss_d    # filtered mechanical power [W]
        e_gov(t)                   # governor error [W]
        x_gov(t)  = 0.0            # governor lag state
        x_int(t)  = τ_o_ss_d      # gate integrator state
        tau_o(t)                   # gate opening [-] (algebraic)
    end
    @equations begin
        # ── Francis affinity turbine (FA.1–FA.3) ──────────────────────────
        Q            ~  tau_o * K_q * D_runner^2 * sqrt(abs(H) + 1e-6)
        eta          ~  eta_max * (1.0 - c_eta * (Q / Q_rated - 1.0)^2)
        P_mech       ~  rho * g * Q * H * eta
        tau_turbine  ~  P_mech / (abs(omega) + 1e-6)

        # ── Generator: load torque (P_load is STEPPED by callback) ──────
        # GA.1 with D_d synchronising torque
        tau_gen      ~  (P_load / omega_ref) *
                        (1.0 + D_d_gen * (omega - omega_ref) / omega_ref)

        # ── Rotor swing equation (M.1) ──────────────────────────────────
        D(omega)     ~  (tau_turbine - tau_gen) / J

        # ── GGOV1 droop governor (GG.1–GG.5) ────────────────────────────
        # Eq. GG.1 – load reference with physical-unit droop
        #   P_ref = P_set + K_droop · (ω_ref − ω)
        #   K_droop = P_r/(σ·ω_ref)  ensures ΔP [W] ↔ Δω [rad/s]
        P_ref        ~  P_set + K_droop * (omega_ref - omega)

        # Eq. GG.2 – first-order power measurement lag
        D(P_meas)    ~  (P_mech - P_meas) / T_p

        # Eq. GG.3 – governor error
        e_gov        ~  P_ref - P_meas

        # Eq. GG.4 – governor lag (first-order, time const = T_gov)
        D(x_gov)     ~  (e_gov - x_gov) / T_gov

        # Gate integrator
        D(x_int)     ~  K_i_gov * e_gov

        # Eq. GG.5 – gate position (clamped to [0, 1])
        tau_o        ~  clamp(x_gov + x_int, 0.0, 1.0)
    end
end

# ── Compile and set up ODE ────────────────────────────────────────────────────
@mtkcompile sys = TrollheimDroop()

u0 = [
    sys.omega   => ω_ref_d,
    sys.P_meas  => P_mech_ss_d,
    sys.x_gov   => 0.0,
    sys.x_int   => τ_o_ss_d,
]

# P_load is stepped at t=50; P_set (governor setpoint) stays constant
p0 = [sys.P_load => P_mech_ss_d]

tspan = (0.0, 100.0)
prob  = ODEProblem(sys, merge(Dict(u0), Dict(p0)), tspan)

# ── 10 % load step at t = 50 s ────────────────────────────────────────────────
# Only the generator LOAD increases; the governor setpoint P_set is unchanged
step_cb = DiscreteCallback(
    (u, t, integ) -> t == 50.0,
    integ -> (integ.ps[sys.P_load] = 1.1 * P_mech_ss_d; nothing),
    save_positions = (true, true)
)

println("Solving Droop simulation (100 s)...")
@time sol = solve(prob, Rodas5P(),
                  tstops   = [50.0],
                  callback = step_cb,
                  abstol   = 1e-7,
                  reltol   = 1e-7)

println("  Retcode  : ", sol.retcode)
@printf "  ω(t=0)   = %.5f rad/s\n" sol[sys.omega][1]
@printf "  ω(t=49)  = %.5f rad/s\n" sol(49.0; idxs=sys.omega)
@printf "  ω(t=55)  = %.5f rad/s\n" sol(55.0; idxs=sys.omega)
@printf "  ω(t=100) = %.5f rad/s\n" sol[sys.omega][end]
@printf "  Δω_ss    = %.5f rad/s  (predicted %.5f rad/s)\n"  sol[sys.omega][end]-ω_ref_d Δω_predicted
@printf "  Δω%%      = %.3f %%\n" (sol[sys.omega][end]-ω_ref_d)/ω_ref_d*100
println()

# ── Result arrays ─────────────────────────────────────────────────────────────
tv      = sol.t
ω_v     = sol[sys.omega]
τo_v    = sol[sys.tau_o]
Pm_v    = sol[sys.P_mech] ./ 1e6
Δω_v    = (ω_v .- ω_ref_d) ./ ω_ref_d .* 100.0   # % deviation

# ── Plot ──────────────────────────────────────────────────────────────────────
default(size=(900, 700), linewidth=2,
        tickfontsize=10, guidefontsize=11, legendfontsize=9, margin=5Plots.mm)

p1 = plot(tv, ω_v,
          label="ω(t)  [rad/s]", ylabel="Speed  [rad/s]",
          legend=:topleft, color=:steelblue)
hline!([ω_ref_d],                    ls=:dash, lc=:red,    lw=1.5, label="ω_ref")
hline!([ω_ref_d + Δω_predicted],     ls=:dot,  lc=:orange, lw=1.5,
       label="ω_ss predicted  (-1 %)")
vline!([50.0], ls=:dot, lc=:black, lw=1, label="Load step  t=50 s")

p2 = plot(tv, τo_v,
          label="τ_o(t)", ylabel="Gate opening  [pu]",
          ylims=(0.7, 1.05), legend=:topleft, color=:darkorange)
vline!([50.0], ls=:dot, lc=:black, lw=1, label=false)

p3 = plot(tv, Pm_v,
          label="P_mech  [MW]", ylabel="Power  [MW]",
          legend=:topleft, color=:seagreen)
hline!([P_mech_ss_d/1e6], ls=:dash, lc=:grey, lw=1.5,
       label="P0 = $(round(P_mech_ss_d/1e6,digits=1)) MW")
hline!([1.1*P_mech_ss_d/1e6], ls=:dash, lc=:red, lw=1.5,
       label="P_step = $(round(1.1*P_mech_ss_d/1e6,digits=1)) MW")
vline!([50.0], ls=:dot, lc=:black, lw=1, label=false)

p4 = plot(tv, Δω_v,
          label="Δω/ω_ref  [%]", ylabel="Freq. deviation  [%]",
          legend=:topleft, color=:firebrick)
hline!([0.0], ls=:dash, lc=:black, lw=1, label="Zero deviation")
hline!([Δω_predicted/ω_ref_d*100], ls=:dot, lc=:orange, lw=1.5,
       label="Predicted -1 %")
vline!([50.0], ls=:dot, lc=:black, lw=1, label="t = 50 s")

ptop = plot(p1, p2, p3, p4,
            layout = (2, 2),
            plot_title = "Trollheim HPP — Droop Governor  σ=0.1  (10 % Load Step at t=50 s)",
            xlabel = "Time  [s]",
            link = :x)

savefig(ptop, joinpath(@__DIR__, "images", "droop_trollheim_step.png"))
println("Plot saved  →  trollheim_droop_control/images/droop_trollheim_step.png")
