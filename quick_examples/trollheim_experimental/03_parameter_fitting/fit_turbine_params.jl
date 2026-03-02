# =============================================================================
# fit_turbine_params.jl  —  Strategy B: Steady-State Parameter Fitting
# Trollheim HPP
#
# Method:
#   At the measured steady-state operating window (t=1200-2400s), fit the
#   turbine model parameters [K_q, eta_max, c_eta] analytically so that
#   the model best reproduces measured (H, Q, P_mech, eta) over all SS points.
#
# Model equations being fitted (FrancisTurbineAffinity):
#   Q     = tau_o * K_q * D_runner^2 * sqrt(H)        (FA.1)
#   eta   = eta_max * (1 - c_eta*(Q/Q_rated - 1)^2)   (FA.2)
#   P_mech= rho * g * Q * H * eta                      (FA.3)
#
# Fitting strategy (no Optim needed — purely analytical least squares):
#
#   K_q:             OLS on FA.1 — single-parameter linear regression
#                    K_q = sum(f_i * Q_i) / sum(f_i^2)
#                    where  f_i = tau_o_i * D_runner^2 * sqrt(H_i)
#
#   [eta_max, c_eta]: Linearise FA.2 as:
#                    eta_i = eta_max - eta_max*c_eta * x_i^2
#                         = a + b*(-x_i^2)
#                    where  x_i = Q_i/Q_rated - 1
#                    Then OLS: [ones  -x^2] \ eta_vec  -> [a, b]
#                    eta_max = a,  c_eta = b/a
#
# Inputs:
#   ../01_preprocessing/trollheim_processed.csv
#
# Outputs:
#   images/parameter_fit.png
#   fitted_params.toml   (plug these back into openloop_validation.jl / sim scripts)
# =============================================================================

using CSV, DataFrames, Statistics
using Plots
using Printf

# ── 1. Load data ───────────────────────────────────────────────────────────────
csv_path = joinpath(@__DIR__, "..", "01_preprocessing", "trollheim_processed.csv")
df       = CSV.read(csv_path, DataFrame)

# ── 2. Physical constants ─────────────────────────────────────────────────────
const rho      = 1000.0
const g_acc    = 9.81
const D_runner = 2.5
const ETA_GEN  = 0.985            # must match preprocess_trollheim.jl

# ── 3. Extract steady-state window ────────────────────────────────────────────
const T_SS_START = 1200.0
const T_SS_END   = 2400.0

ss = filter(row -> T_SS_START <= row.t <= T_SS_END && !isnan(row.eta), df)

Q_ss     = Float64.(ss.Q)
H_ss     = Float64.(ss.H)
tau_o_ss = Float64.(ss.tau_o)
eta_ss   = Float64.(ss.eta)         # measured: P_mech/(rho*g*Q*H)
P_mech_ss= Float64.(ss.P_mech)

n_pts = length(Q_ss)
println("=== Strategy B: Steady-State Parameter Fitting ===")
@printf "  Window: t = %.0f s -> %.0f s  (%d points)\n" T_SS_START T_SS_END n_pts
println()

# ── 4. Fit K_q (from FA.1) ────────────────────────────────────────────────────
# Q_i = K_q * f_i   where   f_i = tau_o_i * D_runner^2 * sqrt(H_i)
f_vec   = tau_o_ss .* D_runner^2 .* sqrt.(H_ss)       # [m3/s] / K_q
K_q_fit = sum(f_vec .* Q_ss) / sum(f_vec .^ 2)        # OLS 1-param

Q_model_ss = K_q_fit .* f_vec                          # model-predicted Q [m3/s]

# ── 5. Fit [eta_max, c_eta] (from FA.2) ───────────────────────────────────────
# Use Q_rated = mean measured SS flow (efficiency peak at SS operating point)
Q_rated_fit = mean(Q_ss)

x_vec = Q_model_ss ./ Q_rated_fit .- 1.0              # (Q/Q_rated - 1), should ≈ 0 at SS
A_mat = hcat(ones(n_pts), -(x_vec .^ 2))              # [1  -(Q/Qr-1)^2] matrix
coefs = A_mat \ eta_ss                                 # OLS: [eta_max;  eta_max*c_eta]

eta_max_fit = coefs[1]
c_eta_fit   = coefs[2] / eta_max_fit

eta_model_ss = eta_max_fit .* (1.0 .- c_eta_fit .* (x_vec .^ 2))
P_mech_model_ss = rho .* g_acc .* Q_model_ss .* H_ss .* eta_model_ss

# ── 6. Fit quality metrics ────────────────────────────────────────────────────
rmse_Q   = sqrt(mean((Q_model_ss     .- Q_ss).^2))
rmse_eta = sqrt(mean((eta_model_ss   .- eta_ss).^2))
rmse_P   = sqrt(mean((P_mech_model_ss .- P_mech_ss).^2)) / 1e6

println("Fitted parameters:")
println("  K_q       = $(round(K_q_fit,      sigdigits = 6))")
println("  eta_max   = $(round(eta_max_fit,   sigdigits = 6))")
println("  c_eta     = $(round(c_eta_fit,     sigdigits = 6))")
println("  Q_rated   = $(round(Q_rated_fit,   sigdigits = 5))  m3/s")
println()
println("Fit quality (SS window):")
@printf "  RMSE Q      = %.4f m3/s  (%.2f %%)\n"  rmse_Q  (rmse_Q/mean(Q_ss)*100)
@printf "  RMSE eta    = %.5f       (%.2f %%)\n"  rmse_eta (rmse_eta/mean(eta_ss)*100)
@printf "  RMSE P_mech = %.3f MW    (%.2f %%)\n"  rmse_P   (rmse_P/mean(P_mech_ss)*1e6*100)
println()
println("Mean SS operating point:")
@printf "  H      = %.1f m\n"    mean(H_ss)
@printf "  Q      = %.2f m3/s\n" mean(Q_ss)
@printf "  tau_o  = %.4f pu\n"   mean(tau_o_ss)
@printf "  P_mech = %.1f MW\n"   mean(P_mech_ss)/1e6
@printf "  eta    = %.4f\n"      mean(eta_ss)

# ── 7. Save fitted_params.toml ────────────────────────────────────────────────
toml_path = joinpath(@__DIR__, "fitted_params.toml")
open(toml_path, "w") do io
    println(io, "# Trollheim HPP — fitted turbine parameters")
    println(io, "# Source: steady-state window t=$(Int(T_SS_START))–$(Int(T_SS_END)) s")
    println(io, "# Generated by: fit_turbine_params.jl (Strategy B)")
    println(io, "")
    println(io, "[turbine]")
    println(io, "K_q       = $(round(K_q_fit,      sigdigits = 6))")
    println(io, "eta_max   = $(round(eta_max_fit,   sigdigits = 6))")
    println(io, "c_eta     = $(round(c_eta_fit,     sigdigits = 6))")
    println(io, "Q_rated   = $(round(Q_rated_fit,   sigdigits = 5))   # m3/s")
    println(io, "D_runner  = $(D_runner)              # m  (fixed from turbine geometry)")
    println(io, "")
    println(io, "[fit_quality]")
    println(io, "rmse_Q_m3s   = $(round(rmse_Q,   sigdigits = 4))")
    println(io, "rmse_eta     = $(round(rmse_eta, sigdigits = 4))")
    println(io, "rmse_P_MW    = $(round(rmse_P,   sigdigits = 4))")
    println(io, "")
    println(io, "[measured_operating_point]")
    println(io, "H_mean_m     = $(round(mean(H_ss),          digits = 2))")
    println(io, "Q_mean_m3s   = $(round(mean(Q_ss),          digits = 3))")
    println(io, "tau_o_mean   = $(round(mean(tau_o_ss),      digits = 4))")
    println(io, "P_elec_MW    = $(round(mean(P_mech_ss)*ETA_GEN/1e6, digits = 2))")
    println(io, "eta_hydraulic= $(round(mean(eta_ss),        digits = 4))")
end
println("Saved -> 03_parameter_fitting/fitted_params.toml")

# ── 8. Plots ──────────────────────────────────────────────────────────────────
t_ss = Float64.(ss.t)

# Panel 1: Q fit over time
p_q = plot(t_ss, Q_ss;
           label  = "Measured Q",
           color  = :black, lw = 1.5, ls = :dash,
           ylabel = "Q  [m3/s]",
           title  = "Strategy B — Parameter Fit Quality  |  Trollheim HPP")
plot!(p_q, t_ss, Q_model_ss;
      label = "FA.1 fit  (K_q=$(round(K_q_fit, sigdigits=4)))",
      color = :steelblue, lw = 1.5)

# Panel 2: efficiency hill curve
q_norm_range = range(0.70, 1.10, 200)
eta_curve    = eta_max_fit .* (1.0 .- c_eta_fit .* (q_norm_range .- 1.0).^2)

p_eta = scatter(Q_ss ./ Q_rated_fit, eta_ss;
                label  = "Measured eta",
                ms = 2, ma = 0.25, color = :black,
                xlabel = "Q / Q_rated",
                ylabel = "Hydraulic efficiency",
                title  = "FA.2 efficiency hill curve fit")
plot!(p_eta, q_norm_range, eta_curve;
      label = "Fitted: eta_max=$(round(eta_max_fit, sigdigits=4)),  c_eta=$(round(c_eta_fit, sigdigits=3))",
      color = :steelblue, lw = 2.0)

# Panel 3: P_mech fit over time
p_p = plot(t_ss, P_mech_ss ./ 1e6;
           label  = "Measured P_mech",
           color  = :black, lw = 1.5, ls = :dash,
           ylabel = "P_mech  [MW]")
plot!(p_p, t_ss, P_mech_model_ss ./ 1e6;
      label = "FA.3 model",
      color = :steelblue, lw = 1.5)

ptop = plot(p_q, p_eta, p_p;
            layout        = (3, 1),
            size          = (900, 900),
            legend        = :topleft,
            left_margin   = 8Plots.mm,
            bottom_margin = 3Plots.mm)

mkpath(joinpath(@__DIR__, "images"))
out_png = joinpath(@__DIR__, "images", "parameter_fit.png")
savefig(ptop, out_png)
println("Plot saved -> 03_parameter_fitting/images/parameter_fit.png")
