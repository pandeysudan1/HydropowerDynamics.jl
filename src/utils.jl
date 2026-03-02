# ============================================================================ #
# utils.jl  –  System-level equations and helper functions
#
# Equations from Section 6 & 7 of the HydroMTK Mathematical Reference.
# ============================================================================ #

# ── Fluid constants (defaults) ───────────────────────────────────────────────
const RHO_0  = 1000.0    # kg/m³  reference density
const B_BULK = 2.0e9     # Pa     bulk modulus of water
const MU_W   = 1.0e-3    # Pa·s   dynamic viscosity
const P_ATM  = 101_325.0 # Pa     atmospheric pressure
const G_GRAV = 9.81      # m/s²   gravitational acceleration

# ── 6.1  Gross and Net Head  [m] ─────────────────────────────────────────────
"""
    gross_head(Z_reservoir, Z_tailwater)

Return H_gross = Z_reservoir - Z_tailwater  [m].  (Eq. SL.1)
"""
gross_head(Z_res, Z_tail) = Z_res - Z_tail

"""
    net_head(H_gross, h_f_penstock, h_f_draft = 0.0)

Return H_net = H_gross - h_f_penstock - h_f_draft  [m].  (Eq. SL.2)
"""
net_head(H_gross, h_f_p, h_f_d = 0.0) = H_gross - h_f_p - h_f_d

"""
    darcy_head_loss(f_D, L, D, v, g = G_GRAV)

Darcy-Weisbach friction head loss:  h_f = f_D*(L/D)*v²/(2g)  [m].  (Eq. SL.3)
"""
darcy_head_loss(f_D, L, D, v, g = G_GRAV) = f_D * (L / D) * v^2 / (2g)

# ── 6.2  Power Cascade  [W] ──────────────────────────────────────────────────
"""
    power_cascade(rho, g, Q, H_gross, H_net, eta_t, eta_gen, eta_xfmr)

Return a `NamedTuple` with the full power cascade:
  `P_hydraulic`, `P_available`, `P_turbine`, `P_electrical`, `P_grid`
and overall `eta_plant`.  (Eqs. PC.1–PC.7)
"""
function power_cascade(rho, g, Q, H_gross, H_net, eta_t, eta_gen, eta_xfmr = 1.0)
    P_hydraulic  = rho * g * Q * H_gross
    P_available  = rho * g * Q * H_net
    P_turbine    = P_available * eta_t
    P_electrical = P_turbine  * eta_gen
    P_grid       = P_electrical * eta_xfmr
    eta_plant    = P_grid / P_hydraulic
    return (; P_hydraulic, P_available, P_turbine, P_electrical, P_grid, eta_plant)
end

"""
    plant_efficiency(eta_hyd, eta_mech, eta_gen, eta_xfmr = 1.0)

Overall plant efficiency as product of component efficiencies.  (Eq. PC.7)
"""
plant_efficiency(eta_hyd, eta_mech, eta_gen, eta_xfmr = 1.0) =
    eta_hyd * eta_mech * eta_gen * eta_xfmr

# ── 6.3  Water Hammer  ───────────────────────────────────────────────────────
"""
    wave_speed(B, rho, D, E_pipe, e_wall)

Pressure-wave celerity in a buried pipeline:
    a = sqrt(B / (rho * (1 + B*D/(E_pipe*e_wall))))  [m/s].  (Eq. WH.2)
"""
wave_speed(B, rho, D, E_pipe, e_wall) =
    sqrt(B / (rho * (1 + B * D / (E_pipe * e_wall))))

"""
    joukowsky_pressure(rho, a_wave, Delta_v)

Maximum instantaneous pressure rise from gate closure:
    Δp = ρ · a · Δv  [Pa].  (Eq. WH.1)
"""
joukowsky_pressure(rho, a_wave, Delta_v) = rho * a_wave * abs(Delta_v)

"""
    critical_closure_time(L, a_wave)

Critical gate-closure time  T_c = 2L/a  [s].  (Eq. WH.3)
Closures faster than T_c experience full Joukowsky pressure.
"""
critical_closure_time(L, a_wave) = 2L / a_wave

# ── 3.1  Dimensionless turbine parameters ────────────────────────────────────
"""
    unit_speed(n_rpm, D, H)

Unit speed  n₁₁ = n·D/√H  [rpm·m^(1/2)].  (Eq. T.1)
"""
unit_speed(n, D, H) = n * D / sqrt(H)

"""
    unit_discharge(Q, D, H)

Unit discharge  Q₁₁ = Q/(D²·√H)  [(m^(1/2))/s].  (Eq. T.2)
"""
unit_discharge(Q, D, H) = Q / (D^2 * sqrt(H))

"""
    hydraulic_efficiency(P_mech, rho, g, Q, H)

η = P_mech / (ρ·g·Q·H).  (Eq. T.4)
"""
hydraulic_efficiency(P_mech, rho, g, Q, H) = P_mech / (rho * g * Q * H)
