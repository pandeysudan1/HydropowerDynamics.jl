# ============================================================================ #
# mechanical.jl  –  Rotating mass and generator models
#
# Section 4 of the HydroMTK Mathematical Reference.
#
# Components
#   RotorInertia    – lumped J model for turbine + shaft + generator rotor
#   SimpleGenerator – algebraic torque-speed characteristic (Section 4.3)
# ============================================================================ #

# ============================================================================ #
# RotorInertia  (Section 4.1)
#
#   J * dω/dt = τ_turbine - τ_generator - τ_friction   (M.1)
#   dθ/dt     = ω                                       (M.2)
#   τ_friction = b_v * ω
# ============================================================================ #
"""
    RotorInertia(; name, J, b_v, omega_0)

Lumped rotating mass (turbine runner + shaft + generator rotor).

**Parameters**
- `J`       : polar moment of inertia [kg·m²]   (default 5000.0)
- `b_v`     : viscous damping [N·m·s/rad]        (default 5.0)
- `omega_0` : initial angular velocity [rad/s]   (default 0.0)

**Ports**
- `flange_turbine`   : RotationalPort – torque input from turbine
- `flange_generator` : RotationalPort – torque output to generator
"""
@mtkmodel RotorInertia begin
    @parameters begin
        J       = 5000.0,  [description = "Polar inertia [kg·m²]"]
        b_v     = 5.0,     [description = "Viscous damping [N·m·s/rad]"]
        omega_0 = 0.0,     [description = "Initial angular velocity [rad/s]"]
    end
    @variables begin
        omega(t) = 0.0,  [description = "Angular velocity [rad/s]"]
        theta(t) = 0.0,  [description = "Angular position [rad]"]
    end
    @components begin
        flange_turbine   = RotationalPort()   # +τ in  (torque driver)
        flange_generator = RotationalPort()   # +τ out (load)
    end
    @equations begin
        # Both flanges rotate at the same speed
        flange_turbine.omega   ~ omega
        flange_generator.omega ~ omega
        flange_turbine.phi     ~ theta
        flange_generator.phi   ~ theta

        # Eq. M.2 – kinematics
        D(theta) ~ omega

        # Eq. M.1 – angular momentum (friction = b_v * ω)
        J * D(omega) ~
            flange_turbine.tau  -           # turbine drives rotor
            (-flange_generator.tau) -       # generator reacts (flow sign flip)
            b_v * omega                     # bearing friction
    end
end

# ============================================================================ #
# SimpleGenerator  (Section 4.3)
#
# Algebraic torque-speed characteristic:
#   τ_gen = (P_rated / ω_rated) * (1 + D_d*(ω - ω_s)/ω_s)   (GA.1)
#   P_elec = τ_gen * ω * η_gen                                (GA.2)
# ============================================================================ #
"""
    SimpleGenerator(; name, P_rated, omega_rated, omega_s, D_d, eta_gen)

Simplified algebraic generator – suitable for long-term stability and
control design studies.

**Parameters**
- `P_rated`      : rated electrical power [W]         (default 50e6)
- `omega_rated`  : rated angular velocity [rad/s]     (default 157.08 ≈ 1500 rpm)
- `omega_s`      : synchronous angular velocity [rad/s] (default 157.08)
- `D_d`          : damping coefficient [-]             (default 1.5)
- `eta_gen`      : generator efficiency [-]            (default 0.97)

**Port**    : `flange` (RotationalPort – input from rotor)
**Outputs** : `tau_gen` [N·m],  `P_elec` [W]
"""
@mtkmodel SimpleGenerator begin
    @parameters begin
        P_rated     = 50.0e6,  [description = "Rated power [W]"]
        omega_rated = 157.08,  [description = "Rated speed [rad/s]"]
        omega_s     = 157.08,  [description = "Synchronous speed [rad/s]"]
        D_d         = 1.5,     [description = "Damping coefficient"]
        eta_gen     = 0.97,    [description = "Generator efficiency"]
    end
    @variables begin
        tau_gen(t),   [description = "Generator reaction torque [N·m]"]
        P_elec(t),    [description = "Electrical output power [W]"]
    end
    @components begin
        flange = RotationalPort()
    end
    @equations begin
        # Eq. GA.1 – algebraic torque-speed
        tau_gen ~ (P_rated / omega_rated) *
                  (1.0 + D_d * (flange.omega - omega_s) / omega_s)

        # Eq. GA.2 – electrical power
        P_elec  ~ tau_gen * abs(flange.omega) * eta_gen

        # Convention: generator reaction torque opposes rotation
        flange.tau ~ -tau_gen
    end
end
