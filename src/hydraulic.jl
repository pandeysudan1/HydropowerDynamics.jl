# ============================================================================ #
# hydraulic.jl  –  Hydraulic component models
#
# All equations reference the HydroMTK Mathematical Reference, Section 2.
#
# Components
#   Reservoir   – pressure head boundary condition
#   Penstock    – 1-D water-hammer pipe with inertia and friction
#   SurgeTank   – free-surface pressure buffer
#   DraftTube   – kinetic energy recovery duct
#   GuideVane   – variable-opening orifice (flow control)
# ============================================================================ #

# ────────────────────────────────────────────────────────────────────────────
# 2.1  Fluid properties (defaults; overridden per component via @parameters)
# ────────────────────────────────────────────────────────────────────────────
#   rho_0  = 1000 kg/m³   reference density
#   B      = 2e9  Pa       bulk modulus
#   mu     = 1e-3 Pa·s     dynamic viscosity
#   p_atm  = 101 325 Pa    atmospheric pressure
#   g      = 9.81 m/s²     gravity

# ============================================================================ #
# Reservoir  (Section 2.2)
#
# Equations:
#   p_port = p_atm + rho * g * H          (Res.1)
#   dm     enforced by Kirchhoff at the node (Res.2)
# ============================================================================ #
"""
    Reservoir(; name, H, rho, g, p_atm)

Upstream reservoir providing a constant pressure head boundary.

**Parameters**
- `H`     : water head above port centreline [m]   (default 100.0)
- `rho`   : water density [kg/m³]                  (default 1000.0)
- `g`     : gravity [m/s²]                         (default 9.81)
- `p_atm` : atmospheric pressure [Pa]              (default 101 325)

**Port** : `port` (HydraulicPort – outlet)
"""
@mtkmodel Reservoir begin
    @parameters begin
        H    = 100.0,   [description = "Head above port [m]"]
        rho  = 1000.0,  [description = "Water density [kg/m³]"]
        g    = 9.81,    [description = "Gravity [m/s²]"]
        p_atm = 101_325.0, [description = "Atm. pressure [Pa]"]
    end
    @components begin
        port = HydraulicPort()
    end
    @equations begin
        # Eq. Res.1 – pressure at port is hydrostatic head
        port.p ~ p_atm + rho * g * H
    end
end

# ============================================================================ #
# Penstock  (Section 2.3)
#
# State variables: dm  (mass flow [kg/s]),  p_avg (average pipe pressure [Pa])
#
# Momentum eq (P.1):
#   L/A * d(dm)/dt = p_in - p_out - Δp_f
#
# Darcy-Weisbach friction (P.2, P.4, P.5):
#   Δp_f = f_D*(L/D)*(dm/(ρ·A))²·ρ/2
#   f_D  = 0.25 / (log10(e/(3.7D) + 5.74/Re^0.9))^2   (Swamee-Jain)
#   Re   = |dm|*D / (μ*A)
#
# Continuity / compressibility (P.3):
#   d(p_avg)/dt = B/(ρ·Vol) * (dm_in - dm_out)
# ============================================================================ #
"""
    Penstock(; name, L, D_pipe, rho, B, mu, e)

1-D penstock with water-hammer inertia and Darcy-Weisbach friction.

**Parameters**
- `L`      : pipe length [m]               (default 500.0)
- `D_pipe` : internal diameter [m]         (default 2.0)
- `rho` : water density [kg/m³]         (default 1000.0)
- `B`   : bulk modulus [Pa]             (default 2.0e9)
- `mu`  : dynamic viscosity [Pa·s]      (default 1e-3)
- `e`   : pipe roughness [m]            (default 1.5e-5)

**Ports** : `port_a` (inlet), `port_b` (outlet)
"""
@mtkmodel Penstock begin
    @parameters begin
        L      = 500.0,   [description = "Pipe length [m]"]
        D_pipe = 2.0,     [description = "Internal diameter [m]"]
        rho    = 1000.0,  [description = "Water density [kg/m³]"]
        B      = 2.0e9,   [description = "Bulk modulus [Pa]"]
        mu     = 1.0e-3,  [description = "Dynamic viscosity [Pa·s]"]
        e      = 1.5e-5,  [description = "Pipe roughness [m]"]
    end
    @variables begin
        dm(t) = 0.0,    [description = "Mass flow rate [kg/s]"]
        p_avg(t) = 0.0, [description = "Average pipe pressure [Pa]"]
        f_D(t),         [description = "Darcy friction factor [-]"]
        Re(t),          [description = "Reynolds number [-]"]
        dp_f(t),        [description = "Friction pressure drop [Pa]"]
    end
    @components begin
        port_a = HydraulicPort()   # inlet
        port_b = HydraulicPort()   # outlet
    end
    @equations begin
        # ── Derived geometry
        # A  = π D_pipe² / 4   (scalar constant embedded below)
        # Vol = A * L

        # ── Reynolds number (Eq. P.5)
        Re ~ abs(dm) * D_pipe / (mu * (π * D_pipe^2 / 4))

        # ── Multi-regime Darcy friction factor (Eq. P.4)  [laminar + transition + turbulent]
        f_D ~ darcy_factor(Re, D_pipe, e)

        # ── Darcy friction pressure drop (Eq. P.2)
        dp_f ~ f_D * (L / D_pipe) * (dm / (rho * (π * D_pipe^2 / 4)))^2 * rho / 2

        # ── Mass / momentum balance at connectors
        port_a.dm + port_b.dm ~ 0.0            # Kirchhoff flow (K2)
        port_a.dm ~ dm                         # sign convention: dm>0 → enters port_a

        # ── Momentum equation – water hammer inertia (Eq. P.1)
        D(dm) ~ (π * D_pipe^2 / 4) / L * (port_a.p - port_b.p - dp_f)

        # ── Compressibility continuity (Eq. P.3)
        # K2 sets port_b.dm = -port_a.dm, so (port_a.dm + port_b.dm) = 0 → p_avg stays at IC
        # Using + here (not -) prevents the sign error that would make D(p_avg) = 2*B*dm/(ρ·A·L)
        D(p_avg) ~ B / (rho * (π * D_pipe^2 / 4) * L) * (port_a.dm + port_b.dm)
    end
end

# ============================================================================ #
# SurgeTank  (Section 2.4)
#
# Free-surface buffer; couples hydraulic port to water level Z.
#
#   A_t * dZ/dt    = dm_port / rho     (ST.1)
#   p_port         = p_atm + rho*g*Z  (ST.2)
# ============================================================================ #
"""
    SurgeTank(; name, A_t, Z_0, D_riser, L_riser, e_riser, rho, g, p_atm)

Surge tank – free-surface pressure attenuator with riser throttle friction.

The riser shaft connecting the tunnel to the tank chamber introduces a
flow-dependent pressure drop that damps surge oscillations. The roughness
height `e_riser` is the primary tuning parameter for damping strength.

**Parameters**
- `A_t`     : tank chamber cross-sectional area [m²]   (default 9.08  — Trollheim: π/4·3.4²)
- `Z_0`     : initial water level [m]                   (default 370.0 — Trollheim SS level)
- `D_riser` : riser shaft diameter [m]                  (default 3.4   — Trollheim geometry)
- `L_riser` : riser shaft length [m]                    (default 87.0  — Trollheim height extent)
- `e_riser` : riser roughness height [m]                (default 1e-2  — unlined rock)
- `rho`     : water density [kg/m³]                     (default 1000.0)
- `g`       : gravity [m/s²]                            (default 9.81)
- `p_atm`   : atmospheric pressure [Pa]                 (default 101 325)

**Port** : `port` (HydraulicPort) — connects at base of riser

**Notes**
- `e_riser = 4.6e-5` m  →  smooth steel  (minimal damping)
- `e_riser = 1.0e-3` m  →  concrete/shotcrete  (moderate damping)
- `e_riser = 1.0e-2` m  →  unlined blasted rock  (strong damping)
"""
@mtkmodel SurgeTank begin
    @parameters begin
        A_t     = 9.08,      [description = "Tank chamber cross-section [m²]"]
        Z_0     = 370.0,     [description = "Initial water level [m]"]
        D_riser = 3.4,       [description = "Riser shaft diameter [m]"]
        L_riser = 87.0,      [description = "Riser shaft length [m]"]
        e_riser = 1.0e-2,    [description = "Riser roughness height [m]"]
        rho     = 1000.0,    [description = "Water density [kg/m³]"]
        g       = 9.81,      [description = "Gravity [m/s²]"]
        p_atm   = 101_325.0, [description = "Atm. pressure [Pa]"]
        mu      = 1.0e-3,    [description = "Dynamic viscosity [Pa·s]"]
    end
    @variables begin
        Z(t) = 370.0,        [description = "Water surface elevation in tank [m]"]
        Re_riser(t),         [description = "Reynolds number in riser [-]"]
        f_riser(t),          [description = "Darcy friction factor in riser [-]"]
        v_riser(t),          [description = "Flow velocity in riser [m/s]"]
        dp_riser(t),         [description = "Riser friction pressure drop [Pa]"]
    end
    @components begin
        port = HydraulicPort()
    end
    @equations begin
        # Riser geometry (area = π/4 * D_riser²)
        v_riser  ~ port.dm / (rho * (π * D_riser^2 / 4.0))

        # Reynolds number in riser
        Re_riser ~ abs(port.dm) * D_riser / (mu * (π * D_riser^2 / 4.0))

        # Multi-regime friction factor (laminar + transition + turbulent)
        f_riser  ~ darcy_factor(Re_riser, D_riser, e_riser)

        # Darcy-Weisbach pressure drop across riser (Eq. ST.3, sign-preserving)
        dp_riser ~ f_riser * (L_riser / D_riser) * rho / 2.0 * v_riser * abs(v_riser)

        # Eq. ST.1 – free-surface continuity
        D(Z) ~ port.dm / (rho * A_t)

        # Eq. ST.2 – hydrostatic pressure at port (base of riser) including throttle loss
        # sign: inflow (dm>0) raises Z → riser friction reduces available pressure at port
        port.p ~ p_atm + rho * g * Z - dp_riser
    end
end

# ============================================================================ #
# DraftTube  (Section 2.5)
#
# Diverging duct with friction and exit loss.
#
#   p_in - p_out = -ρ/2*(v_out² - v_in²) + Δp_f + Δp_exit   (DT.1)
#   Δp_exit = C_exit * ρ/2 * v_out²                          (DT.2)
# ============================================================================ #
"""
    DraftTube(; name, A_in, A_out, L, f_D, C_exit, rho)

Draft tube – kinetic energy recovery after the turbine runner.

**Parameters**
- `A_in`   : inlet cross-section [m²]        (default 0.785)
- `A_out`  : outlet cross-section [m²]       (default 3.14)
- `L`      : tube length [m]                 (default 6.0)
- `f_D`    : Darcy friction factor [-]       (default 0.02)
- `C_exit` : exit loss coefficient [-]       (default 1.0)
- `rho`    : water density [kg/m³]           (default 1000.0)

**Ports** : `port_a` (inlet), `port_b` (outlet)
"""
@mtkmodel DraftTube begin
    @parameters begin
        A_in   = 0.785,   [description = "Inlet area [m²]"]
        A_out  = 3.14,    [description = "Outlet area [m²]"]
        L      = 6.0,     [description = "Tube length [m]"]
        f_D    = 0.02,    [description = "Darcy friction factor"]
        C_exit = 1.0,     [description = "Exit loss coefficient"]
        rho    = 1000.0,  [description = "Water density [kg/m³]"]
    end
    @variables begin
        v_in(t),   [description = "Inlet velocity [m/s]"]
        v_out(t),  [description = "Outlet velocity [m/s]"]
        dp_f(t),   [description = "Friction pressure drop [Pa]"]
    end
    @components begin
        port_a = HydraulicPort()   # turbine exit / draft tube inlet
        port_b = HydraulicPort()   # tailrace outlet
    end
    @equations begin
        # Sign convention: dm > 0 enters port_a, exits port_b
        port_a.dm + port_b.dm ~ 0.0

        # Velocities from mass flow and density
        v_in  ~ port_a.dm / (rho * A_in)
        v_out ~ port_a.dm / (rho * A_out)

        # Darcy friction in the expanding duct (averaged velocity)
        dp_f  ~ f_D * (L / sqrt(A_in)) * ((v_in + v_out) / 2)^2 * rho / 2

        # Eq. DT.1 + DT.2  (energy equation)
        port_a.p - port_b.p ~
            -rho / 2 * (v_out^2 - v_in^2) + dp_f + C_exit * rho / 2 * v_out^2
    end
end

# ============================================================================ #
# GuideVane  (Section 2.6)  –  variable-opening orifice
#
#   dm = C_d * A_0 * tau_o * sqrt(2*rho*|Δp|) * sign(Δp)   (GV.1)
#   A_0 = π/4 * D_0²                                         (GV.2)
# ============================================================================ #
"""
    GuideVane(; name, D_0, C_d, rho)

Variable-opening guide vane / wicket gate orifice.

**Parameters**
- `D_0`  : reference diameter at full opening [m]   (default 1.0)
- `C_d`  : discharge coefficient [-]                 (default 0.7)
- `rho`  : water density [kg/m³]                     (default 1000.0)

**Ports**
- `port_a` : upstream hydraulic port
- `port_b` : downstream hydraulic port
- `opening` : signal input – normalised gate opening τₒ ∈ [0, 1]
"""
@mtkmodel GuideVane begin
    @parameters begin
        D_0 = 1.0,    [description = "Reference diameter at full open [m]"]
        C_d = 0.7,    [description = "Discharge coefficient"]
        rho = 1000.0, [description = "Water density [kg/m³]"]
    end
    @variables begin
        dm(t) = 0.0,    [description = "Mass flow rate [kg/s]"]
        Delta_p(t),     [description = "Pressure drop across vane [Pa]"]
    end
    @components begin
        port_a  = HydraulicPort()
        port_b  = HydraulicPort()
        opening = SignalInPort()
    end
    @equations begin
        port_a.dm + port_b.dm ~ 0.0          # Kirchhoff
        port_a.dm ~ dm

        Delta_p ~ port_a.p - port_b.p

        # Eq. GV.1 + GV.2  (orifice equation, sign-preserving)
        dm ~ C_d * (π / 4 * D_0^2) * opening.u *
             sqrt(2 * rho * abs(Delta_p) + 1e-6) * sign(Delta_p + 1e-20)
    end
end
