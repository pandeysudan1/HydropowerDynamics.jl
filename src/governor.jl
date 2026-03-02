# ============================================================================ #
# governor.jl  –  Speed governor and control models
#
# Section 5 of the HydroMTK Mathematical Reference.
#
# PIDGovernor   – IEEE-style PID speed governor  (Section 5.1)
# GGOV1Governor – IEEE Std 1207 GGOV1 model      (Section 5.2)
# ============================================================================ #

# ============================================================================ #
# PIDGovernor  (Section 5.1)
#
#   e(t)     = ω_ref - ω(t)                              (Gov.1)
#   u_PID(t) = K_p*e + K_i*∫e + K_d*de/dt               (Gov.2)
#   τ_o(t)   = clamp(u_PID, 0, 1)                        (Gov.3)
#   |dτ_o/dt| ≤ R_open  (opening) / R_close (closing)   (Gov.4)
#
# Note: Rate limits are implemented as soft bounds on D(tau_o).
# ============================================================================ #
"""
    PIDGovernor(; name, K_p, K_i, K_d, omega_ref, R_open, R_close, tau_0)

IEEE-style PID speed governor.

**Parameters**
- `K_p`      : proportional gain               (default 2.0)
- `K_i`      : integral gain [1/s]             (default 1.0)
- `K_d`      : derivative gain [s]             (default 0.1)
- `omega_ref`: speed reference [rad/s]         (default 157.08)
- `R_open`   : max gate opening rate [1/s]     (default 0.1)
- `R_close`  : max gate closing rate [1/s]     (default 0.2)
- `tau_0`    : initial gate position [-]        (default 0.8)

**Ports**
- `speed_in`  : SignalInPort  – measured rotor speed ω
- `gate_out`  : SignalOutPort – gate opening command τₒ ∈ [0, 1]
"""
@mtkmodel PIDGovernor begin
    @parameters begin
        K_p       = 2.0,    [description = "Proportional gain"]
        K_i       = 1.0,    [description = "Integral gain [1/s]"]
        K_d       = 0.1,    [description = "Derivative gain [s]"]
        omega_ref = 157.08, [description = "Speed reference [rad/s]"]
        R_open    = 0.1,    [description = "Max gate open rate [1/s]"]
        R_close   = 0.2,    [description = "Max gate close rate [1/s]"]
        tau_0     = 0.8,    [description = "Initial gate opening"]
    end
    @variables begin
        e(t),       [description = "Speed error [rad/s]"]
        e_int(t) = 0.0, [description = "Integral of speed error"]
        de_dt(t),   [description = "Speed error derivative [rad/s²]"]
        u_pid(t),   [description = "PID output (before saturation)"]
        tau_o(t) = 0.8, [description = "Gate opening output [-]"]
    end
    @components begin
        speed_in = SignalInPort()
        gate_out = SignalOutPort()
    end
    @equations begin
        # Eq. Gov.1 – speed error
        e ~ omega_ref - speed_in.u

        # Integral of error
        D(e_int) ~ e

        # Derivative of error
        de_dt ~ D(e)

        # Eq. Gov.2 – PID output
        u_pid ~ K_p * e + K_i * e_int + K_d * de_dt

        # Eq. Gov.3 & Gov.4 – soft rate-limited saturation
        # Gate position is a state; its rate is limited and clamped
        D(tau_o) ~ max(-R_close, min(R_open,
                   (clamp(u_pid, 0.0, 1.0) - tau_o) / 0.01))

        # Drive output port
        gate_out.u ~ tau_o
    end
end

# ============================================================================ #
# GGOV1Governor  (Section 5.2)
#
# IEEE Std 1207 GGOV1 model with droop, power measurement lag, governor lag.
#
#   P_ref  = P_set + (ω_ref - ω) / R_droop    (GG.1)
#   T_p * d(P_meas)/dt = P_mech - P_meas       (GG.2)
#   e_gov  = P_ref - P_meas                    (GG.3)
#   T_gov * d(x_gov)/dt = e_gov - x_gov        (GG.4)
#   τ_o   = clamp(x_gov + x_int, τ_min, τ_max) (GG.5)
# ============================================================================ #
"""
    GGOV1Governor(; name, R_droop, T_p, T_gov, omega_ref, P_set,
                    tau_min, tau_max, K_i, tau_0)

IEEE Std 1207 GGOV1 governor with droop, power lag and gate integrator.

**Parameters**
- `R_droop`  : speed droop (regulation) [-]           (default 0.05)
- `T_p`      : power measurement time constant [s]    (default 0.05)
- `T_gov`    : governor lag time constant [s]         (default 0.4)
- `omega_ref`: speed reference [rad/s]                (default 157.08)
- `P_set`    : load setpoint [pu or W]                (default 1.0)
- `tau_min`  : minimum gate opening [-]               (default 0.0)
- `tau_max`  : maximum gate opening [-]               (default 1.0)
- `K_i`      : gate integrator gain [1/s]             (default 1.0)
- `tau_0`    : initial gate opening [-]               (default 0.8)

**Ports**
- `speed_in`  : SignalInPort  – rotor speed ω [rad/s]
- `power_in`  : SignalInPort  – measured mechanical power P_mech [W or pu]
- `gate_out`  : SignalOutPort – gate opening τₒ ∈ [τ_min, τ_max]
"""
@mtkmodel GGOV1Governor begin
    @parameters begin
        R_droop  = 0.05,   [description = "Speed droop"]
        T_p      = 0.05,   [description = "Power measurement lag [s]"]
        T_gov    = 0.4,    [description = "Governor lag [s]"]
        omega_ref = 157.08,[description = "Speed reference [rad/s]"]
        P_set    = 1.0,    [description = "Load setpoint [pu or W]"]
        tau_min  = 0.0,    [description = "Min gate opening"]
        tau_max  = 1.0,    [description = "Max gate opening"]
        K_i      = 1.0,    [description = "Gate integrator gain [1/s]"]
        tau_0    = 0.8,    [description = "Initial gate opening"]
    end
    @variables begin
        P_ref(t),          [description = "Load reference with droop [pu/W]"]
        P_meas(t) = 0.0,   [description = "Filtered mechanical power [pu/W]"]
        e_gov(t),          [description = "Governor error"]
        x_gov(t) = 0.0,    [description = "Governor lag state"]
        x_int(t) = 0.8,    [description = "Gate integrator state"]
        tau_o(t) = 0.8,    [description = "Gate opening [-]"]
    end
    @components begin
        speed_in = SignalInPort()
        power_in = SignalInPort()
        gate_out = SignalOutPort()
    end
    @equations begin
        # Eq. GG.1 – load reference with droop
        P_ref ~ P_set + (omega_ref - speed_in.u) / R_droop

        # Eq. GG.2 – power measurement first-order lag
        D(P_meas) ~ (power_in.u - P_meas) / T_p

        # Eq. GG.3 – governor error
        e_gov ~ P_ref - P_meas

        # Eq. GG.4 – governor lag (first-order)
        D(x_gov) ~ (e_gov - x_gov) / T_gov

        # Gate integrator
        D(x_int) ~ K_i * e_gov

        # Eq. GG.5 – gate position (clamped)
        tau_o ~ clamp(x_gov + x_int, tau_min, tau_max)

        # Output
        gate_out.u ~ tau_o
    end
end
