# ============================================================================ #
# governor.jl  –  Speed governor and control models
#
# Section 5 of the HydroMTK Mathematical Reference.
#
# Built on ModelingToolkitStandardLibrary.Blocks (MTSL) as the causal layer.
# Dynamic components (lags, integrators, rate-limiters, PID) come directly
# from MTSL so that analysis-point tagging and get_sensitivity() work out of
# the box.  Custom @equations handle only static algebraic terms (droop sum,
# error signal) that do not have MTSL equivalents with parameter-driven inputs.
#
# MTSL blocks used:
#   Blocks.Constant          – constant signal source (ω_ref setpoint)
#   Blocks.FirstOrder        – first-order lag  k/(sT+1)
#   Blocks.Integrator        – pure integrator  k/s
#   Blocks.Derivative        – filtered derivative  ks/(sT+1)
#   Blocks.Add               – signal summer with gain weights
#   Blocks.Limiter           – static saturation  [y_min, y_max]
#   Blocks.SlewRateLimiter   – rate-limited output  (dτ ≤ R_open / R_close)
#   Blocks.LimPID            – complete PID with anti-windup and output limits
#
# PIDGovernor  (Section 5.1) – ISO-style PID speed governor
# GGOV1Governor(Section 5.2) – IEEE Std 1207 GGOV1 droop governor
# ============================================================================ #

# ============================================================================ #
# PIDGovernor  (Section 5.1)
#
# Implements the IEEE ISO PID speed governor using MTSL's LimPID block:
#
#   e(t)     = ω_ref - ω(t)                                    (Gov.1)
#   u_PID(t) = K_p*(e + ∫e/T_i + T_d·de/dt)                  (Gov.2)
#   τ_o      = sat(u_PID, 0, 1)           via LimPID u_max/min (Gov.3)
#   dτ_o/dt ∈ [-R_close, R_open]          via SlewRateLimiter  (Gov.4)
#
# MTSL block topology:
#
#   Constant(ω_ref) ──┐
#                     ├─► LimPID ──► SlewRateLimiter ──► gate_out
#   speed_in ─────────┘  (reference)  (measurement)
#
# LimPID wires:  reference = Constant(ω_ref),  measurement = speed_in
# ============================================================================ #
"""
    PIDGovernor(; name, K_p, T_i, T_d, omega_ref, R_open, R_close, tau_0)

ISO-style PID speed governor built from MTSL `Blocks.LimPID` + `Blocks.SlewRateLimiter`.

**Parameters**
- `K_p`      : proportional gain                    (default 2.0)
- `T_i`      : integral time constant [s]   (= 1/Kᵢ) (default 1.0)
- `T_d`      : derivative time constant [s] (= Kd/Kp)(default 0.1)
- `omega_ref`: speed reference [rad/s]              (default 157.08)
- `R_open`   : max gate opening rate [pu/s]          (default 0.1)
- `R_close`  : max gate closing rate [pu/s]          (default 0.2)
- `tau_0`    : initial gate position [-]             (default 0.8)

**Ports** (MTSL `Blocks.RealInput` / `RealOutput`)
- `speed_in`  : RealInput  – measured rotor speed ω [rad/s]
- `gate_out`  : RealOutput – gate opening command τₒ ∈ [0, 1]

**Internal MTSL blocks** (accessible as `gov.pid.*`, `gov.slew.*`)
- `ref`  : `Blocks.Constant`          – ω_ref set-point
- `pid`  : `Blocks.LimPID`            – PID + anti-windup + output limits
- `slew` : `Blocks.SlewRateLimiter`   – gate rate limiter

**Backward-compat alias**: `tau_o` variable mirrors `slew.output.u`
"""
@mtkmodel PIDGovernor begin
    @parameters begin
        K_p       = 2.0,    [description = "Proportional gain"]
        T_i       = 1.0,    [description = "Integral time constant [s]  (= 1/K_i)"]
        T_d       = 0.1,    [description = "Derivative time constant [s] (= K_d/K_p)"]
        omega_ref = 157.08, [description = "Speed reference [rad/s]"]
        R_open    = 0.1,    [description = "Max gate opening rate [pu/s]"]
        R_close   = 0.2,    [description = "Max gate closing rate [pu/s]"]
        tau_0     = 0.8,    [description = "Initial gate position [-]"]
    end
    @variables begin
        tau_o(t) = 0.8,  [description = "Gate opening output (alias of slew.output.u) [-]"]
    end
    @components begin
        # ── External ports (MTSL RealInput / RealOutput) ─────────────────────
        speed_in = Blocks.RealInput()
        gate_out = Blocks.RealOutput()

        # ── Constant speed setpoint ───────────────────────────────────────────
        # Provides ω_ref as a signal so LimPID.reference can be connected.
        ref  = Blocks.Constant(k = omega_ref)

        # ── LimPID – PID controller with output saturation + anti-windup ─────
        # reference = ω_ref (set-point),  measurement = ω (plant output)
        # u_max / u_min clamp the raw PID output to [0, 1] before slew limiting.
        # Nd sets derivative filter bandwidth (reasonable: 8–20).
        pid  = Blocks.LimPID(
                   k     = K_p,
                   Ti    = T_i,
                   Td    = T_d,
                   Nd    = 10.0,
                   u_max = 1.0,
                   u_min = 0.0)

        # ── SlewRateLimiter – gate rate limits (Gov.4) ───────────────────────
        # rising  =  R_open  [pu/s]   (max opening speed)
        # falling = -R_close [pu/s]   (max closing speed, must be ≤ 0)
        # Td is the internal derivative filter time constant (keep small).
        slew = Blocks.SlewRateLimiter(
                   rising  = R_open,
                   falling = -R_close,
                   Td      = 0.001,
                   y_start = tau_0)
    end
    @equations begin
        # ── Connect LimPID: ref → reference, speed_in → measurement ──────────
        connect(ref.output,     pid.reference)
        connect(speed_in,       pid.measurement)

        # ── Connect rate limiter ──────────────────────────────────────────────
        connect(pid.ctr_output, slew.input)

        # ── Drive output port ─────────────────────────────────────────────────
        connect(slew.output,    gate_out)

        # ── Backward-compat alias for simulation scripts ──────────────────────
        tau_o ~ gate_out.u
    end
end

# ============================================================================ #
# GGOV1Governor  (Section 5.2)
#
# IEEE Std 1207 GGOV1 model rebuilt using MTSL dynamic blocks.
#
# Equation layout (unchanged from math reference):
#   P_ref  = P_set + (ω_ref - ω) / R_droop          (GG.1, algebraic)
#   T_p · dP_meas/dt = P_mech - P_meas              (GG.2) ← FirstOrder pow_lag
#   e_gov  = P_ref - pow_lag.output.u               (GG.3, algebraic)
#   T_gov · dx_gov/dt = e_gov - x_gov               (GG.4) ← FirstOrder gov_lag
#   dx_int/dt = K_i · e_gov                          (gate integrator) ← Integrator
#   tau_o  = clamp(x_gov + x_int, τ_min, τ_max)     (GG.5) ← Add + Limiter
#
# MTSL block topology:
#
#   power_in ──► [FirstOrder T_p] ──► pow_lag ─────────────────────┐
#                                                                    ↓
#   (P_ref = P_set + (ω_ref-ω)/R_droop) ──► e_gov (algebraic) ──► [Add]
#                                                                    │
#                                             ┌─── [FirstOrder T_gov] ─→ x_gov ─┐
#                                             │                                   ↓
#                                             └─── [Integrator K_i]  ─→ x_int ──► [Add] → [Limiter] → gate_out
#
# GG.1 and GG.3 remain @equations (static, parameter-driven, no MTSL equivalent
# for a Constant that uses a tunable parameter P_set inside @components).
# ============================================================================ #
"""
    GGOV1Governor(; name, R_droop, T_p, T_gov, omega_ref, P_set,
                    tau_min, tau_max, K_i, tau_0)

IEEE Std 1207 GGOV1 governor rebuilt with MTSL blocks for every dynamic element.

**Parameters**
- `R_droop`  : speed droop [rad·s/W]       (= σ·ω_ref/P_r)  (default 0.05)
- `T_p`      : power measurement lag [s]                      (default 0.05)
- `T_gov`    : governor lag [s]                                (default 0.4)
- `omega_ref`: speed reference [rad/s]                        (default 157.08)
- `P_set`    : load setpoint [W]            (tunable)          (default 1.0)
- `tau_min`  : minimum gate opening [-]                        (default 0.0)
- `tau_max`  : maximum gate opening [-]                        (default 1.0)
- `K_i`      : gate integrator gain [1/s]                      (default 1.0)
- `tau_0`    : initial gate opening [-]                        (default 0.8)

**Ports** (MTSL `Blocks.RealInput` / `RealOutput`)
- `speed_in`  : RealInput  – rotor speed ω [rad/s]
- `power_in`  : RealInput  – turbine mechanical power P_mech [W]
- `gate_out`  : RealOutput – gate opening τₒ ∈ [τ_min, τ_max]

**Internal MTSL blocks** (accessible as `gov.pow_lag.*` etc.)
- `pow_lag`  : `Blocks.FirstOrder`  – GG.2 power measurement filter
- `gov_lag`  : `Blocks.FirstOrder`  – GG.4 governor response lag
- `gate_int` : `Blocks.Integrator`  – gate position integrator
- `gate_sum` : `Blocks.Add`         – GG.5 sum  x_gov + x_int
- `gate_lim` : `Blocks.Limiter`     – GG.5 clamp to [τ_min, τ_max]

**Backward-compat aliases**: `P_meas` mirrors `pow_lag.output.u`;
`x_gov` mirrors `gov_lag.output.u`; `x_int` mirrors `gate_int.output.u`;
`tau_o` mirrors `gate_lim.output.u`.
"""
@mtkmodel GGOV1Governor begin
    @parameters begin
        R_droop   = 0.05,   [description = "Speed droop [rad·s/W]"]
        T_p       = 0.05,   [description = "Power measurement lag [s]"]
        T_gov     = 0.4,    [description = "Governor lag [s]"]
        omega_ref = 157.08, [description = "Speed reference [rad/s]"]
        P_set     = 1.0,    [description = "Load setpoint [W]", tunable = true]
        tau_min   = 0.0,    [description = "Min gate opening"]
        tau_max   = 1.0,    [description = "Max gate opening"]
        K_i       = 1.0,    [description = "Gate integrator gain [1/s]"]
        tau_0     = 0.8,    [description = "Initial gate opening"]
    end
    @variables begin
        # ── Static algebraic intermediates ──────────────────────────────────
        P_ref(t),               [description = "Power reference with droop [W]"]
        e_gov(t),               [description = "Governor error signal [W]"]
        # ── Backward-compatibility aliases (mirror internal MTSL block states) ─
        P_meas(t) = 0.0,        [description = "Filtered P_mech (= pow_lag.output.u) [W]"]
        x_gov(t)  = 0.0,        [description = "Governor lag output (= gov_lag.output.u)"]
        x_int(t)  = 0.8,        [description = "Gate integrator state (= gate_int.output.u)"]
        tau_o(t)  = 0.8,        [description = "Gate opening (= gate_lim.output.u) [-]"]
    end
    @components begin
        # ── External ports ───────────────────────────────────────────────────
        speed_in = Blocks.RealInput()
        power_in = Blocks.RealInput()
        gate_out = Blocks.RealOutput()

        # ── GG.2 Power measurement first-order lag  T_p·dP_meas/dt = P_mech - P_meas
        pow_lag  = Blocks.FirstOrder(T = T_p,   k = 1.0)

        # ── GG.4 Governor first-order lag  T_gov·dx_gov/dt = e_gov - x_gov
        gov_lag  = Blocks.FirstOrder(T = T_gov, k = 1.0)

        # ── Gate integrator  dx_int/dt = K_i · e_gov
        gate_int = Blocks.Integrator(k = K_i)

        # ── GG.5a Gate position sum  s = x_gov + x_int
        gate_sum = Blocks.Add()

        # ── GG.5b Gate position clamp  τ_o = clamp(s, τ_min, τ_max)
        gate_lim = Blocks.Limiter(y_max = tau_max, y_min = tau_min)
    end
    @equations begin
        # ── GG.1: load reference with speed droop (algebraic) ─────────────────
        # P_set is a tunable parameter; ω_ref is a fixed parameter.
        # Both are constants in the signal sense — no MTSL block needed.
        P_ref ~ P_set + (omega_ref - speed_in.u) / R_droop

        # ── GG.2: power measurement lag ──────────────────────────────────────
        connect(power_in, pow_lag.input)

        # ── GG.3: governor error (algebraic) ─────────────────────────────────
        e_gov ~ P_ref - pow_lag.output.u

        # ── GG.4: governor lag input driven by e_gov (algebraic injection) ───
        gov_lag.input.u  ~ e_gov

        # ── Gate integrator input driven by e_gov ────────────────────────────
        gate_int.input.u ~ e_gov

        # ── GG.5: gate computation chain ─────────────────────────────────────
        connect(gov_lag.output,  gate_sum.input1)
        connect(gate_int.output, gate_sum.input2)
        connect(gate_sum.output, gate_lim.input)
        connect(gate_lim.output, gate_out)

        # ── Backward-compatibility aliases ────────────────────────────────────
        P_meas ~ pow_lag.output.u
        x_gov  ~ gov_lag.output.u
        x_int  ~ gate_int.output.u
        tau_o  ~ gate_lim.output.u
    end
end
