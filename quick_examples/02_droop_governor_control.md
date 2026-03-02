# Quick Example 2 — Droop Governor Control
## Trollheim HPP (Operating Point B)

**Plant:** Trollheim Hydropower Plant (high-head Francis turbine)  
**Operating point:** $H_n = 371\ \text{m}$, $\dot{V}_n = 37\ \text{m}^3/\text{s}$, $P_r \approx 130\ \text{MW}$, $n = 375\ \text{RPM}$  
**Scenario:** 10 % load increase at $t = 50\ \text{s}$; plant in steady state for $t < 50\ \text{s}$  
**Total simulation:** 100 s

**Simulation script:** [`droop_trollheim_test.jl`](trollheim_droop_control/droop_trollheim_test.jl)

### Step response plot

![Droop step response](trollheim_droop_control/images/droop_trollheim_step.png)

---

## 1. Plant Parameters

| Symbol | Description | Value |
|--------|-------------|-------|
| $H_n$ | Nominal head | 371 m |
| $\dot{V}_n$ | Nominal discharge | 37 m³/s |
| $P_r$ | Rated mechanical power | ≈130 MW |
| $n$ | Nominal speed | 375 RPM |
| $L_p$ | Penstock length | 500 m |
| $A_p$ | Penstock cross-section | $4\pi$ m² |
| $T_w$ | Water time constant | 0.44 s |
| $T_{gs}$ | Gate servo time constant | 0.2 s |
| $T_{ps}$ | Pilot servo time constant | 1.75 s |
| $T_d$ | Dashpot valve time constant | 0.2 s |
| $\sigma$ | Static droop coefficient | 0.1 |
| $\delta$ | Transient droop coefficient | 0.04 |
| $\dot{u}_\text{max}$ | Max gate opening rate | 0.05 pu/s |
| $\dot{u}_\text{min}$ | Max gate closing rate | 0.2 pu/s |
| $J$ | Rotor inertia | $7 \times 10^5$ kg·m² |

---

## 2. Droop Governor Theory

A **droop governor** allows a proportional steady-state speed deviation in response to load changes. This enables stable **load sharing** among multiple units connected to the same grid.

### 2.1 Speed Droop Definition

The static droop $\sigma$ (regulation) relates the steady-state speed change to the change in gate opening:

$$\sigma = \frac{\Delta \omega / \omega_\text{ref}}{\Delta \tau_o} \quad \Rightarrow \quad \Delta \tau_o = \frac{\Delta \omega}{\sigma \cdot \omega_\text{ref}}$$

With $\sigma = 0.1$ (from Trollheim benchmark), a **1 % frequency drop** causes a **10 % gate opening increase**.

### 2.2 GGOV1 Load Reference with Droop

The GGOV1 model (IEEE Std 1207) computes a load reference that incorporates speed droop:

$$P_\text{ref}(t) = P_\text{set} + \frac{\omega_\text{ref} - \omega(t)}{R_\text{droop}}$$

where $R_\text{droop} = \sigma \cdot \omega_\text{ref} / P_r$ (in per-unit power base). In the implementation:

$$R_\text{droop} = \sigma = 0.1$$

### 2.3 Power Measurement Lag

The measured mechanical power is filtered through a first-order lag (prevents noise sensitivity):

$$T_p \frac{dP_\text{meas}}{dt} = P_\text{mech} - P_\text{meas}, \qquad T_p = 0.05\ \text{s}$$

### 2.4 Governor Dynamics

Governor error and lag:

$$e_\text{gov}(t) = P_\text{ref}(t) - P_\text{meas}(t)$$

$$T_\text{gov} \frac{dx_\text{gov}}{dt} = e_\text{gov} - x_\text{gov}, \qquad T_\text{gov} \approx T_{ps} = 1.75\ \text{s}$$

Gate integrator:

$$\frac{dx_\text{int}}{dt} = K_i \cdot e_\text{gov}$$

Gate position (clamped):

$$\tau_o = \text{clamp}(x_\text{gov} + x_\text{int},\ \tau_\text{min},\ \tau_\text{max})$$

### 2.5 Steady-State Speed Deviation

Unlike isochronous control, droop introduces a **permanent** speed offset proportional to load change $\Delta P$:

$$\Delta \omega_\text{ss} = -\sigma \cdot \omega_\text{ref} \cdot \frac{\Delta P}{P_r}$$

For a 10 % load step ($\Delta P = 0.1 \cdot P_r$) with $\sigma = 0.1$:

$$\Delta \omega_\text{ss} = -0.1 \times 39.27 \times 0.1 \approx -0.39\ \text{rad/s} \quad (\approx -1\%\ \text{frequency drop})$$

### 2.6 Turbine Equations (same as Example 1)

$$Q = \tau_o \cdot K_q \cdot D^2 \cdot \sqrt{H}, \qquad K_q \approx 0.095\ \text{m}^{-1/2}\text{s}^{-1}$$

$$\eta = \eta_\text{max} \left[1 - c_\eta \left(\frac{Q}{Q_\text{rated}} - 1\right)^2\right]$$

$$P_\text{mech} = \rho g Q H \eta$$

### 2.7 Rotor Swing Equation

$$J \frac{d\omega}{dt} = \tau_\text{turbine} - \tau_\text{generator} - b_v \omega$$

With mechanical starting time $M_a = J \omega_\text{ref}^2 / P_r \approx 8\ \text{s}$.

---

## 3. Julia Code

```julia
using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots

# ── Trollheim HPP — Operating Point A ────────────────────────────────────────
const H_n       = 460.0
const Q_n       = 24.0
const P_r       = 130e6
const n_rpm     = 375.0
const omega_ref = 2π * n_rpm / 60   # ≈ 39.27 rad/s

# Penstock (benchmark table)
const L_p   = 500.0
const D_p   = 2 * sqrt(π)           # A_p = 4π → D = 2√π ≈ 3.545 m

# Turbine sizing
const D_runner = 3.43
const K_q      = Q_n / (D_runner^2 * sqrt(H_n))

# Rotor: mechanical starting time M_a = 8 s
const J_rotor  = 8.0 * P_r / omega_ref^2

# Droop and governor parameters (Trollheim benchmark)
const sigma     = 0.1      # static droop
const T_ps      = 1.75     # pilot servo / governor lag [s]
const T_gs      = 0.2      # gate servo [s]

@mtkmodel TrollheimDroop begin
    @components begin
        # ── Hydraulic circuit ─────────────────────────────────────────────
        reservoir  = Reservoir(H = H_n)
        penstock   = Penstock(L = L_p, D = D_p)
        turbine    = FrancisTurbineAffinity(
                         D       = D_runner,
                         K_q     = K_q,
                         eta_max = 0.92,
                         c_eta   = 0.25,
                         Q_rated = Q_n,
                     )
        tailwater  = Reservoir(H = 0.0)

        # ── Mechanical train ──────────────────────────────────────────────
        rotor      = RotorInertia(
                         J       = J_rotor,
                         b_v     = 10.0,
                         omega_0 = omega_ref,
                     )
        generator  = SimpleGenerator(
                         P_rated     = P_r,
                         omega_rated = omega_ref,
                         omega_s     = omega_ref,
                         D_d         = 2.0,         # damping coefficient
                         eta_gen     = 0.98,
                     )

        # ── GGOV1 droop governor ──────────────────────────────────────────
        # R_droop = sigma = 0.1 (Trollheim static droop)
        # T_gov   = T_ps  = 1.75 s (pilot servo)
        # P_set   = 1.0 pu (rated load)
        governor   = GGOV1Governor(
                         R_droop   = sigma,
                         T_p       = 0.05,
                         T_gov     = T_ps,
                         omega_ref = omega_ref,
                         P_set     = P_r,
                         tau_min   = 0.0,
                         tau_max   = 1.0,
                         K_i       = 1.0 / T_gs,    # 1/T_gs per gate servo
                         tau_0     = 0.8,
                     )
    end
    @equations begin
        # ── Hydraulic connections ─────────────────────────────────────────
        connect(reservoir.port,   penstock.port_a)
        connect(penstock.port_b,  turbine.port_in)
        connect(turbine.port_out, tailwater.port)

        # ── Mechanical connections ────────────────────────────────────────
        connect(turbine.shaft,          rotor.flange_turbine)
        connect(rotor.flange_generator, generator.flange)

        # ── Droop governor loop ───────────────────────────────────────────
        # Speed and power feeds → governor → gate command
        connect(rotor.speed_out,    governor.speed_in)
        connect(turbine.power_out,  governor.power_in)
        connect(governor.gate_out,  turbine.opening)
    end
end

@named sys = TrollheimDroop()
sys = structural_simplify(sys)

# Initial conditions at 80 % gate, rated speed
u0 = [
    sys.rotor.omega       => omega_ref,
    sys.rotor.theta       => 0.0,
    sys.governor.P_meas   => P_r * 0.8,
    sys.governor.x_gov    => 0.0,
    sys.governor.x_int    => 0.8,
    sys.penstock.dm       => 1000.0 * Q_n,
]

# ── 10 % load step at t = 5 s ────────────────────────────────────────────────
# Simulate by increasing generator load setpoint
function load_step_affect!(integrator)
    integrator[sys.generator.D_d] *= 1.1    # increase damping = more load demand
end
cb = PresetTimeCallback(5.0, load_step_affect!)

tspan = (0.0, 60.0)
prob  = ODEProblem(sys, u0, tspan)
sol   = solve(prob, Rodas5P(), callback=cb, abstol=1e-6, reltol=1e-6)

# ── Plot results ──────────────────────────────────────────────────────────────
omega_ss_drop = -sigma * omega_ref * 0.1    # theoretical steady-state drop

p1 = plot(sol, idxs = sys.rotor.omega,
          label = "ω [rad/s]", ylabel = "Speed (rad/s)",
          title  = "Trollheim — Droop Control: Rotor Speed")
hline!([omega_ref],                      ls=:dash,  lc=:red,   label="ω_ref")
hline!([omega_ref + omega_ss_drop],      ls=:dot,   lc=:orange, label="ω_ss (predicted)")

p2 = plot(sol, idxs = sys.governor.tau_o,
          label = "τₒ [-]", ylabel = "Gate opening [-]",
          title  = "Gate Position")

p3 = plot(sol, idxs = sys.turbine.P_mech ./ 1e6,
          label = "P_mech [MW]", ylabel = "Power (MW)",
          title  = "Mechanical Power Output")
hline!([P_r/1e6], ls=:dash, lc=:green, label="P_rated")

p4 = plot(sol, idxs = (sys.rotor.omega .- omega_ref) ./ omega_ref .* 100,
          label = "Δω/ω_ref [%]", ylabel = "Freq. deviation (%)",
          title  = "Frequency Deviation")
hline!([-sigma*10], ls=:dot, lc=:orange, label="Theoretical Δω_ss")

plot(p1, p2, p3, p4, layout=(2,2), size=(1000,700), xlabel="Time (s)")
savefig("quick_examples/trollheim_droop.png")
println("Simulation complete. Plot saved to quick_examples/trollheim_droop.png")
```

---

## 4. Expected Behaviour

After the 10 % load step at $t = 5\ \text{s}$:

| Quantity | Isochronous (Ex. 1) | Droop (This example) |
|---|---|---|
| Steady-state Δω | **0** (restored) | **−0.39 rad/s** (≈ −1%) |
| Gate response | Aggressive (full integral) | Proportional + lagged |
| Grid suitability | Single unit / island | Multi-unit grid |
| Load sharing | Not inherent | Automatic via droop |

### Key observations:

1. **Speed drops** after the load step and reaches a **new steady-state** below $\omega_\text{ref}$ — this is intentional droop behaviour.
2. The **governor lag** $T_\text{gov} = 1.75\ \text{s}$ produces a slower, smoother gate response than the isochronous PID.
3. The **power settles** at the new higher load level without oscillation.
4. Multiple units with the same droop $\sigma = 0.1$ share load **proportionally to their rating**.

### Droop Load Sharing Rule (multi-unit):

For two units with ratings $P_{r,1}$ and $P_{r,2}$ and the same droop $\sigma$:

$$\frac{\Delta P_1}{\Delta P_2} = \frac{P_{r,1}}{P_{r,2}}$$

> **Note:** To restore frequency to $\omega_\text{ref}$ after a droop response, a secondary controller (AGC — Automatic Generation Control) must adjust $P_\text{set}$ over a slower timescale (tens of seconds to minutes). This is outside the scope of this example.
