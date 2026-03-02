# ============================================================================ #
# turbine.jl  –  Turbine models
#
# Section 3 of the HydroMTK Mathematical Reference.
#
# FrancisTurbineAffinity  – quadratic affinity-law model     (Section 3.3)
# PeltonTurbine           – jet-impact model                 (Section 3.4)
# ============================================================================ #

# ============================================================================ #
# FrancisTurbineAffinity  (Section 3.3)
#
# Simplified quadratic characteristic for linearization and control design:
#
#   Q     = τ_o * K_q * D² * √H                              (FA.1)
#   η     = η_max * [1 - c_η*(Q/Q_rated - 1)²]               (FA.2)
#   P_mech= ρ * g * Q * H * η                                 (FA.3)
#   τ_shaft = P_mech / ω                                      (from T.4)
#   dm    = ρ * Q                                             (from F.4)
# ============================================================================ #
"""
    FrancisTurbineAffinity(; name, D, K_q, eta_max, c_eta, Q_rated, rho, g)

Francis turbine – affinity-law (quadratic) performance model.
Suitable for control design and linearization.

**Parameters**
- `D`       : runner diameter [m]            (default 2.5)
- `K_q`     : flow coefficient at rated [-]  (default 0.35)
- `eta_max` : peak efficiency [-]            (default 0.92)
- `c_eta`   : parabolic loss coefficient [-] (default 0.25)
- `Q_rated` : rated volumetric flow [m³/s]   (default 100.0)
- `rho`     : water density [kg/m³]          (default 1000.0)
- `g`       : gravity [m/s²]                 (default 9.81)

**Ports**
- `port_in`  : HydraulicPort – high-pressure inlet
- `port_out` : HydraulicPort – low-pressure outlet
- `opening`  : SignalInPort  – normalised gate opening τₒ ∈ [0, 1]
- `shaft`    : RotationalPort – mechanical output
"""
@mtkmodel FrancisTurbineAffinity begin
    @parameters begin
        D       = 2.5,    [description = "Runner diameter [m]"]
        K_q     = 0.35,   [description = "Flow coefficient at rated"]
        eta_max = 0.92,   [description = "Peak efficiency"]
        c_eta   = 0.25,   [description = "Parabolic loss coefficient"]
        Q_rated = 100.0,  [description = "Rated volumetric flow [m³/s]"]
        rho     = 1000.0, [description = "Water density [kg/m³]"]
        g       = 9.81,   [description = "Gravity [m/s²]"]
    end
    @variables begin
        H(t),         [description = "Net head across turbine [m]"]
        Q(t),         [description = "Volumetric flow rate [m³/s]"]
        eta(t),       [description = "Hydraulic efficiency [-]"]
        P_mech(t),    [description = "Mechanical power output [W]"]
        tau_shaft(t), [description = "Shaft torque [N·m]"]
        dm(t),        [description = "Mass flow rate [kg/s]"]
    end
    @components begin
        port_in  = HydraulicPort()
        port_out = HydraulicPort()
        opening  = SignalInPort()
        shaft    = RotationalPort()
    end
    @equations begin
        # Mass balance at hydraulic ports
        port_in.dm + port_out.dm ~ 0.0
        port_in.dm ~ dm

        # Net head (Eq. T.4 prerequisite)
        H ~ (port_in.p - port_out.p) / (rho * g)

        # Eq. FA.1 – flow from gate opening and head
        Q ~ opening.u * K_q * D^2 * sqrt(abs(H) + 1e-6)

        # Mass flow
        dm ~ rho * Q

        # Eq. FA.2 – parabolic efficiency hill
        eta ~ eta_max * (1.0 - c_eta * (Q / Q_rated - 1.0)^2)

        # Eq. FA.3 – mechanical power
        P_mech ~ rho * g * Q * abs(H) * eta

        # Shaft torque (from Eq. T.4)
        tau_shaft ~ P_mech / (abs(shaft.omega) + 1e-6)

        # Rotational connector: deliver torque, angular velocity from rotor
        shaft.tau ~ tau_shaft
    end
end

# ============================================================================ #
# PeltonTurbine  (Section 3.4)
#
# Jet-impact model:
#
#   v_jet  = C_v * √(2*g*H_net)                                (PL.1)
#   Q      = A_jet(τ_o) * v_jet                                (PL.2)
#   τ_shaft= ρ*Q*v_jet*u₁*(1 + φ*cos β₂) / ω                 (PL.3)
#   u₁     = ω * R_mean                                        (PL.4)
# ============================================================================ #
"""
    PeltonTurbine(; name, R_mean, D_jet_max, C_v, phi, beta2_deg, rho, g)

Pelton turbine – jet-impact performance model.

**Parameters**
- `R_mean`     : pitch circle radius [m]         (default 1.2)
- `D_jet_max`  : max jet diameter (fully open) [m] (default 0.15)
- `C_v`        : velocity coefficient [-]          (default 0.98)
- `phi`        : relative velocity ratio at exit [-] (default 0.85)
- `beta2_deg`  : bucket exit angle [°]             (default 170.0)
- `rho`        : water density [kg/m³]              (default 1000.0)
- `g`          : gravity [m/s²]                     (default 9.81)

**Ports**
- `port_in`  : HydraulicPort – nozzle inlet (high pressure)
- `port_out` : HydraulicPort – bucket outlet (tailrace, low pressure)
- `opening`  : SignalInPort  – needle opening τₒ ∈ [0, 1]
- `shaft`    : RotationalPort – mechanical output
"""
@mtkmodel PeltonTurbine begin
    @parameters begin
        R_mean    = 1.2,    [description = "Pitch circle radius [m]"]
        D_jet_max = 0.15,   [description = "Max jet diameter [m]"]
        C_v       = 0.98,   [description = "Velocity coefficient"]
        phi       = 0.85,   [description = "Relative velocity ratio"]
        beta2_deg = 170.0,  [description = "Bucket exit angle [°]"]
        rho       = 1000.0, [description = "Water density [kg/m³]"]
        g         = 9.81,   [description = "Gravity [m/s²]"]
    end
    @variables begin
        H_net(t),     [description = "Net head at nozzle [m]"]
        v_jet(t),     [description = "Jet velocity [m/s]"]
        A_jet(t),     [description = "Jet area [m²]"]
        Q(t),         [description = "Volumetric flow [m³/s]"]
        dm(t),        [description = "Mass flow rate [kg/s]"]
        u1(t),        [description = "Bucket peripheral velocity [m/s]"]
        tau_shaft(t), [description = "Shaft torque [N·m]"]
    end
    @components begin
        port_in  = HydraulicPort()
        port_out = HydraulicPort()
        opening  = SignalInPort()
        shaft    = RotationalPort()
    end
    @equations begin
        # Mass balance
        port_in.dm + port_out.dm ~ 0.0
        port_in.dm ~ dm

        # Net head (pressure difference converted to head)
        H_net ~ (port_in.p - port_out.p) / (rho * g)

        # Eq. PL.1 – jet velocity
        v_jet ~ C_v * sqrt(2.0 * g * abs(H_net) + 1e-6)

        # Eq. PL.2 – jet area proportional to needle opening²
        #   A_jet = τ_o² * A_max  (area scales as square of opening fraction)
        A_jet ~ opening.u^2 * (π / 4 * D_jet_max^2)

        # Volume and mass flow
        Q  ~ A_jet * v_jet
        dm ~ rho * Q

        # Eq. PL.4 – bucket peripheral velocity
        u1 ~ abs(shaft.omega) * R_mean

        # Eq. PL.3 – shaft torque from jet-bucket impulse
        tau_shaft ~ rho * Q * v_jet * u1 *
                    (1.0 + phi * cosd(beta2_deg)) / (abs(shaft.omega) + 1e-6)

        # Deliver torque to shaft port
        shaft.tau ~ tau_shaft
    end
end
