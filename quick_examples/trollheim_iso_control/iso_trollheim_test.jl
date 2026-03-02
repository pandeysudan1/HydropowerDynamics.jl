# =============================================================================
# iso_trollheim_test.jl
# Trollheim HPP — Isochronous (ISO) PID Governor Step Response
#
# Operating Point B:  H = 371 m,  Q = 37 m³/s,  P ≈ 130 MW,  n = 375 RPM
# Scenario:           10 % load increase applied at t = 50 s
# Simulation:         0 s → 100 s   (steady state before t = 50 s)
#
# Governor type:      Isochronous PID  (no droop, δ = 0)
# Expected outcome:   speed returns exactly to ω_ref after transient
#
# Benchmark data (Trollheim, governor table):
#   T_w   = 0.44 s    R_open  = 0.05 pu/s    R_close = 0.20 pu/s
#   T_gs  = 0.20 s    M_a (mechanical starting time) ≈ 8 s
# =============================================================================

using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots
using Printf

mkpath(joinpath(@__DIR__, "quick_examples", "images"))

# ── Plant constants ───────────────────────────────────────────────────────────
const H_n    = 371.0              # Net head [m]
const Q_n    = 37.0               # Nominal discharge [m³/s]
const n_rpm  = 375.0              # Nominal speed [RPM]
const ω_ref  = 2π * n_rpm / 60   # ≈ 39.27 rad/s
const ρ      = 1000.0
const g_acc  = 9.81

# Turbine sizing
const D_r    = 2.5                # Runner diameter [m]
const η_max_v = 0.97              # Peak efficiency at BEP
const c_eta_v = 0.25              # Parabolic efficiency loss coefficient
const Q_rated_v = Q_n             # BEP at nominal flow (τ_o = 0.90)

# K_q sizing: Q_n = τ_o_0 * K_q * D_r² * √H_n  →  rated gate τ_o_0 = 0.90
const τ_o_ss = 0.90
const K_q_v  = Q_n / (τ_o_ss * D_r^2 * sqrt(H_n))

# Steady-state mechanical power (exact at τ_o = τ_o_ss)
const P_mech_ss = ρ * g_acc * Q_n * H_n * η_max_v   # ≈ 130.6 MW

# Rotor inertia: mechanical starting time M_a = 8 s → J = M_a·P/(ω²)
const J_v    = 8.0 * P_mech_ss / ω_ref^2             # ≈ 6.74×10⁵ kg·m²

# Governor ISO-PID gains (tuned for T_w = 0.44 s)
const Kp = 2.0
const Ki = 2.0
const R_open_v  = 0.05            # max gate opening rate [pu/s]
const R_close_v = 0.20            # max gate closing rate [pu/s]

# Steady-state e_int so that PID output = τ_o_ss at zero speed error
const e_int_ss = τ_o_ss / Ki     # = 0.45

println("=== Trollheim ISO Governor ===")
@printf "  ω_ref     = %.4f rad/s  (%g RPM)\n" ω_ref n_rpm
@printf "  K_q       = %.4f\n" K_q_v
@printf "  P_mech_ss = %.2f MW\n" P_mech_ss/1e6
@printf "  J_rotor   = %.0f kg·m²\n" J_v
println()

# ── Flat MTK model ────────────────────────────────────────────────────────────
@mtkmodel TrollheimISO begin
    @parameters begin
        # Turbine/hydraulic
        H      = H_n
        rho    = ρ
        g      = g_acc
        D_runner = D_r
        K_q    = K_q_v
        eta_max = η_max_v
        c_eta  = c_eta_v
        Q_rated = Q_rated_v
        # Rotor
        J      = J_v
        D_d    = 2.0          # generator synchronising torque coefficient
        # Generator load (stepped at t=50)
        P_rated = P_mech_ss   # [W] – tunable: stepped 10% up by callback
        # ISO-PID governor
        K_p    = Kp
        K_i    = Ki
        omega_ref = ω_ref
        R_open  = R_open_v
        R_close = R_close_v
    end
    @variables begin
        omega(t)  = ω_ref          # rotor angular velocity [rad/s]
        Q(t)                       # volumetric flow [m³/s]
        eta(t)                     # hydraulic efficiency [-]
        P_mech(t)                  # mechanical power [W]
        tau_turbine(t)             # turbine shaft torque [N·m]
        tau_gen(t)                 # generator reaction torque [N·m]
        e(t)                       # speed error [rad/s]
        e_int(t) = e_int_ss        # integrated speed error
        tau_o(t) = τ_o_ss          # gate opening [-]
        u_pid(t)                   # raw PID output
    end
    @equations begin
        # ── Francis affinity turbine (FA.1–FA.3) ──────────────────────────
        Q         ~  tau_o * K_q * D_runner^2 * sqrt(abs(H) + 1e-6)
        eta       ~  eta_max * (1.0 - c_eta * (Q / Q_rated - 1.0)^2)
        P_mech    ~  rho * g * Q * H * eta
        tau_turbine ~ P_mech / (abs(omega) + 1e-6)

        # ── Generator (algebraic, GA.1) ─────────────────────────────────
        tau_gen   ~  (P_rated / omega_ref) *
                     (1.0 + D_d * (omega - omega_ref) / omega_ref)

        # ── Rotor swing equation (M.1) ──────────────────────────────────
        D(omega)  ~  (tau_turbine - tau_gen) / J

        # ── ISO PID governor (Gov.1–Gov.4) ──────────────────────────────
        e         ~  omega_ref - omega
        D(e_int)  ~  e
        u_pid     ~  K_p * e + K_i * e_int
        D(tau_o)  ~  max(-R_close,
                         min(R_open,
                             (clamp(u_pid, 0.0, 1.0) - tau_o) / 0.01))
    end
end

# ── Compile and set up ODE ────────────────────────────────────────────────────
@mtkcompile sys = TrollheimISO()

u0 = [
    sys.omega  => ω_ref,
    sys.e_int  => e_int_ss,
    sys.tau_o  => τ_o_ss,
]

p0 = [sys.P_rated => P_mech_ss]

tspan = (0.0, 100.0)
prob  = ODEProblem(sys, merge(Dict(u0), Dict(p0)), tspan)

# ── 10 % load step at t = 50 s ────────────────────────────────────────────────
step_cb = DiscreteCallback(
    (u, t, integ) -> t == 50.0,
    integ -> (integ.ps[sys.P_rated] = 1.1 * P_mech_ss; nothing),
    save_positions = (true, true)
)

println("Solving ISO simulation (100 s)...")
@time sol = solve(prob, Rodas5P(),
                  tstops    = [50.0],
                  callback  = step_cb,
                  abstol    = 1e-7,
                  reltol    = 1e-7)

println("  Retcode : ", sol.retcode)
@printf "  ω(t=0)  = %.5f rad/s\n" sol[sys.omega][1]
@printf "  ω(t=49) = %.5f rad/s\n" sol(49.0; idxs=sys.omega)
@printf "  ω(t=55) = %.5f rad/s\n" sol(55.0; idxs=sys.omega)
@printf "  ω(t=100)= %.5f rad/s  (vs ref %.5f)\n" sol[sys.omega][end] ω_ref
@printf "  Δω_ss   = %.6f rad/s  (%.4f %%)\n"  sol[sys.omega][end]-ω_ref (sol[sys.omega][end]-ω_ref)/ω_ref*100
println()

# ── Save result arrays ────────────────────────────────────────────────────────
tv      = sol.t
ω_v     = sol[sys.omega]
τo_v    = sol[sys.tau_o]
Pm_v    = sol[sys.P_mech] ./ 1e6
e_v     = ω_v .- ω_ref

# ── Plot ──────────────────────────────────────────────────────────────────────
default(size=(900, 700), linewidth=2,
        tickfontsize=10, guidefontsize=11, legendfontsize=9, margin=5Plots.mm)

p1 = plot(tv, ω_v,
          label="ω(t)  [rad/s]", ylabel="Speed  [rad/s]",
          legend=:topleft, color=:steelblue)
hline!([ω_ref], ls=:dash, lc=:red, lw=1.5, label="ω_ref = $(round(ω_ref,digits=3))")
vline!([50.0], ls=:dot, lc=:black, lw=1, label="Load step  t=50 s")

p2 = plot(tv, τo_v,
          label="τ_o(t)", ylabel="Gate opening  [pu]",
          ylims=(0.7, 1.05), legend=:topleft, color=:darkorange)
vline!([50.0], ls=:dot, lc=:black, lw=1, label=false)

p3 = plot(tv, Pm_v,
          label="P_mech  [MW]", ylabel="Power  [MW]",
          legend=:topleft, color=:seagreen)
hline!([P_mech_ss/1e6], ls=:dash, lc=:grey, lw=1.5,
       label="P0 = $(round(P_mech_ss/1e6,digits=1)) MW")
hline!([1.1*P_mech_ss/1e6], ls=:dash, lc=:red, lw=1.5,
       label="P_step = $(round(1.1*P_mech_ss/1e6,digits=1)) MW")
vline!([50.0], ls=:dot, lc=:black, lw=1, label=false)

p4 = plot(tv, e_v,
          label="e = ω_ref - ω  [rad/s]", ylabel="Speed error  [rad/s]",
          legend=:topleft, color=:firebrick)
hline!([0.0], ls=:dash, lc=:black, lw=1, label="Zero error")
vline!([50.0], ls=:dot, lc=:black, lw=1, label="t = 50 s")

ptop = plot(p1, p2, p3, p4,
            layout = (2, 2),
            plot_title = "Trollheim HPP — Isochronous Governor  (10 % Load Step at t=50 s)",
            xlabel = "Time  [s]",
            link = :x)

savefig(ptop, joinpath(@__DIR__, "images", "iso_trollheim_step.png"))
println("Plot saved  →  trollheim_iso_control/images/iso_trollheim_step.png")
