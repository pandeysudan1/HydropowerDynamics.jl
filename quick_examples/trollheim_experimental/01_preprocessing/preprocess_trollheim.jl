# =============================================================================
# preprocess_trollheim.jl
# Trollheim HPP — Experimental Data Preprocessing & Overview
#
# Source:  quick_examples/data/Trollheim_data.xlsx  (3600 samples, 1 s/sample = 1 hour)
#
# Raw channel -> library variable mapping
# -----------------------------------------------------------------------------
#  Raw channel         Unit   Library var     Description
#  p_penstock_a_  ->  kPa   p_in   [Pa]     Turbine INLET  pressure (spiral casing)
#  p_dt_b_        ->  kPa   p_out  [Pa]     Turbine OUTLET pressure (draft tube)
#  p_surna_a_     ->  moh   h_tailwater [m]  Tailwater head at outlet tunnel
#  tacho_fine_    ->  RPM   omega  [rad/s]   Rotor angular velocity
#  pos_servo_     ->  mm    tau_o  [pu]      Guide vane opening (0=closed, 1=full)
#  power_generator_ -> MW  P_elec [W]       Generator electrical output
#  Q_tunnel_      -> m3/s   Q      [m3/s]   Turbine discharge
#
# Derived variables (matching HydroPowerDynamics.jl library equations)
#   H     = (p_in - p_out) / (rho*g)       [m]    net head across turbine
#           matches: H ~ (port_in.p - port_out.p)/(rho*g)  FrancisTurbineAffinity
#   P_mech= P_elec / eta_gen               [W]    mechanical shaft power
#   eta   = P_mech / (rho*g*Q*H)          [-]    hydraulic efficiency
#           matches: eta ~ eta_max*(1 - c_eta*(Q/Q_rated-1)^2)
#   dm    = rho * Q                        [kg/s]  mass flow (HydraulicPort through var)
#
# Outputs:
#   data/trollheim_processed.csv
#   quick_examples/images/trollheim_overview.png
#   quick_examples/images/trollheim_efficiency.png
# =============================================================================

using XLSX
using DataFrames
using CSV
using Statistics
using Printf
using Plots

mkpath(joinpath(@__DIR__, "..", "images"))

DATA_PATH = joinpath(@__DIR__, "..", "..", "data", "Trollheim_data.xlsx")

# -- Physical constants -------------------------------------------------------
const rho       = 1000.0              # water density      [kg/m3]
const g_acc     = 9.81                # gravity            [m/s2]
const n_nom     = 375.0               # nominal speed      [RPM]
const omega_nom = 2pi * n_nom / 60    # nominal omega      [rad/s]  ~39.27
const GV_MAX_MM = 99.134              # guide vane full-open position [mm]
const ETA_GEN   = 0.985               # assumed generator efficiency [-]

# -- Load raw Excel -----------------------------------------------------------
println("Loading Excel data...")
raw = DataFrame(XLSX.readtable(DATA_PATH, "Data"))
println("  Raw shape: $(size(raw))")
println("  Raw columns: $(names(raw))")

# -- Step 1: rename to library variable names --------------------------------
# HydraulicPort connector   : p [Pa],     dm [kg/s]
# FrancisTurbineAffinity    : H [m],  Q [m3/s], eta [-], P_mech [W], tau_o [pu]
# RotorInertia/RotationalPort: omega [rad/s]
rename!(raw, [
    "p_penstock_a_ Mean "    => "p_in_kPa",      # turbine inlet  pressure [kPa]
    "p_dt_b_ Mean "          => "p_out_kPa",     # turbine outlet pressure [kPa]  (draft tube)
    "p_surna_a_ Mean "       => "h_tailwater_m", # tailwater head [moh = m]
    "tacho_fine_ Mean "      => "n_rpm",          # rotor speed [RPM]
    "pos_servo_ Mean "       => "gv_mm",          # guide vane servo position [mm]
    "power_generator_ Mean " => "P_elec_MW",      # generator electrical power [MW]
    "Q_tunnel_ Mean "        => "Q",              # turbine discharge [m3/s]  (direct match)
])

# -- Step 2: add time axis ---------------------------------------------------
raw[!, :t] = collect(0:(nrow(raw)-1))   # t [s],  1 s resolution
select!(raw, :t, :)                      # move time to first column

# -- Step 3: convert raw units -> library SI units ---------------------------

# Hydraulic pressures: kPa -> Pa
# (HydraulicPort across variable: p [Pa])
raw[!, :p_in]  = raw.p_in_kPa  .* 1e3   # [Pa]  turbine inlet  (spiral casing)
raw[!, :p_out] = raw.p_out_kPa .* 1e3   # [Pa]  turbine outlet (draft tube b)

# Net head across turbine: H = (p_in - p_out) / (rho*g)  [m]
# Equation: H ~ (port_in.p - port_out.p) / (rho * g)  [FrancisTurbineAffinity]
raw[!, :H] = (raw.p_in .- raw.p_out) ./ (rho * g_acc)

# Rotor speed: RPM -> rad/s
# (RotationalPort across variable: omega [rad/s])
raw[!, :omega]    = raw.n_rpm .* (2pi / 60)
raw[!, :omega_pu] = raw.n_rpm ./ n_nom

# Guide vane opening: servo mm -> pu in [0, 1]
# (SignalInPort u variable: tau_o [pu])
raw[!, :tau_o] = clamp.(raw.gv_mm, 0.0, GV_MAX_MM) ./ GV_MAX_MM

# Generator power: MW -> W
raw[!, :P_elec] = raw.P_elec_MW .* 1e6   # [W]

# Mechanical shaft power: P_mech = P_elec / eta_gen  [W]
raw[!, :P_mech] = raw.P_elec ./ ETA_GEN

# Hydraulic efficiency: eta = P_mech / (rho*g*Q*H)  [-]
# Reference model: eta ~ eta_max * (1 - c_eta*(Q/Q_rated - 1)^2)
function calc_eta(P_mech, Q, H)
    P_hydro = rho * g_acc * Q * H   # [W]
    (Q > 0.5 && P_mech > 1e6 && P_hydro > 1e5) ? P_mech / P_hydro : NaN
end
raw[!, :eta] = calc_eta.(raw.P_mech, raw.Q, raw.H)

# Mass flow: dm = rho*Q [kg/s]
# (HydraulicPort through variable: dm [kg/s])
raw[!, :dm] = rho .* raw.Q

# -- Step 4: statistics -------------------------------------------------------
println("\n=== Signal statistics (library variable names) ===")
for (col, lbl, unit) in [
        (:H,         "H",       "m"),
        (:Q,         "Q",       "m3/s"),
        (:omega,     "omega",   "rad/s"),
        (:tau_o,     "tau_o",   "pu"),
        (:P_elec_MW, "P_elec",  "MW"),
        (:eta,       "eta",     "-")]
    v = filter(!isnan, raw[!, col])
    isempty(v) && continue
    println("  $(rpad(lbl,10)) [$(rpad(unit,6))]  min=$(round(minimum(v),digits=3))  mean=$(round(mean(v),digits=3))  max=$(round(maximum(v),digits=3))")
end

# Rated operating point: top-5% by power
P_thresh  = quantile(filter(x -> x > 1.0, raw.P_elec_MW), 0.95)
rated_idx = findall(raw.P_elec_MW .>= P_thresh)
H_r    = mean(raw.H[rated_idx])
Q_r    = mean(raw.Q[rated_idx])
P_r    = mean(raw.P_elec_MW[rated_idx])
n_r    = mean(raw.n_rpm[rated_idx])
tau_r  = mean(raw.tau_o[rated_idx])
eta_rv = filter(!isnan, raw.eta[rated_idx])
eta_r  = isempty(eta_rv) ? NaN : mean(eta_rv)

println("\n=== Rated operating point (top-5% by power) ===")
@printf "  H     = %.1f m\n"                       H_r
@printf "  Q     = %.2f m3/s\n"                    Q_r
@printf "  P_elec= %.2f MW\n"                      P_r
@printf "  omega = %.4f rad/s  (%.2f RPM)\n"       n_r*(2pi/60) n_r
@printf "  tau_o = %.3f pu  (guide vane)\n"        tau_r
@printf "  eta   = %.4f  (hydraulic efficiency)\n" eta_r

# -- Step 5: save processed CSV ----------------------------------------------
out_cols = [:t, :p_in, :p_out, :h_tailwater_m, :H, :Q, :dm,
            :omega, :omega_pu, :tau_o, :P_elec, :P_mech, :eta,
            :n_rpm, :P_elec_MW]
out_csv = joinpath(@__DIR__, "trollheim_processed.csv")
CSV.write(out_csv, raw[!, out_cols])
println("\nProcessed CSV saved -> $out_csv")
println("  Library-aligned columns:")
for c in out_cols
    println("    ", c)
end

# -- Step 6: overview plot (4-panel) -----------------------------------------
tv = raw.t

default(size=(1050, 780), linewidth=1.5,
        tickfontsize=9, guidefontsize=10, legendfontsize=8, margin=5Plots.mm)

p1 = plot(tv, raw.omega,
          ylabel="omega [rad/s]", label="omega(t)", legend=:topleft, color=:steelblue)
hline!([omega_nom], ls=:dash, lc=:red, lw=1.5,
       label="omega_nom = $(round(omega_nom, digits=2)) rad/s")

p2 = plot(tv, raw.tau_o,
          ylabel="tau_o [pu]", label="tau_o(t)", legend=:topleft, color=:darkorange)

p3 = plot(tv, raw.P_elec_MW,
          ylabel="P_elec [MW]", label="P_elec(t)", legend=:topleft, color=:seagreen)
hline!([P_r], ls=:dash, lc=:grey, lw=1.5,
       label="P_rated = $(round(P_r, digits=1)) MW")

p4 = plot(tv, raw.H,
          ylabel="H [m]", label="H = (p_in-p_out)/(rho*g)", legend=:topleft, color=:purple)
hline!([H_r], ls=:dash, lc=:grey, lw=1.5,
       label="H_rated = $(round(H_r, digits=1)) m")

pov = plot(p1, p2, p3, p4,
           layout=(2, 2),
           plot_title="Trollheim HPP  Experimental Data Overview  (3600 s)",
           xlabel="Time [s]", link=:x)

out1 = joinpath(@__DIR__, "..", "images", "trollheim_overview.png")
savefig(pov, out1)
println("\nOverview plot saved -> $out1")

# -- Step 7: efficiency scatter (hill chart proxy) ---------------------------
mask  = .!isnan.(raw.eta) .& (raw.P_elec_MW .> 5.0)
Q_v   = raw.Q[mask]
eta_v = raw.eta[mask]
tau_v = raw.tau_o[mask]

peff = scatter(Q_v, eta_v,
               marker_z=tau_v, color=:turbo, colorbar=true,
               label=false, colorbar_title="tau_o [pu]",
               xlabel="Q [m3/s]",
               ylabel="eta  hydraulic efficiency [-]",
               title="Trollheim HPP  eta vs Q  (colour = tau_o  guide vane opening)",
               legend=:topleft, markersize=3, markerstrokewidth=0,
               size=(820, 520))

out2 = joinpath(@__DIR__, "..", "images", "trollheim_efficiency.png")
savefig(peff, out2)
println("Efficiency scatter saved -> $out2")

println("\nDone.")
