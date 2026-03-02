# Quick Example 1 — Isochronous Governor Control
## Trollheim HPP (Operating Point B)

**Plant:** Trollheim Hydropower Plant (high-head Francis turbine)  
**Operating point:** $H_n = 371\ \text{m}$, $\dot{V}_n = 37\ \text{m}^3/\text{s}$, $P_r \approx 130\ \text{MW}$, $n = 375\ \text{RPM}$  
**Scenario:** 10 % load increase at $t = 50\ \text{s}$; plant in steady state for $t < 50\ \text{s}$  
**Total simulation:** 100 s

**Simulation script:** [`iso_trollheim_test.jl`](trollheim_iso_control/iso_trollheim_test.jl)

### Step response plot

![Isochronous step response](trollheim_iso_control/images/iso_trollheim_step.png)

---

## 1. Plant Parameters

| Symbol | Description | Value |
|--------|-------------|-------|
| $H_n$ | Nominal head | 371 m |
| $\dot{V}_n$ | Nominal discharge | 37 m³/s |
| $P_r$ | Rated mechanical power | ≈130 MW |
| $n$ | Nominal speed | 375 RPM |
| $J$ | Rotor inertia ($M_a = 8\ \text{s}$) | ≈$6.74 \times 10^5$ kg·m² |
| $T_w$ | Water time constant | 0.44 s |
| $T_{gs}$ | Gate servo time constant | 0.2 s |
| $J$ | Rotor inertia | $7 \times 10^5$ kg·m² |

Synchronous angular velocity:

$$\omega_{\text{ref}} = \frac{2\pi \cdot n}{60} = \frac{2\pi \cdot 375}{60} \approx 39.27\ \text{rad/s}$$

Runner diameter (from unit discharge $Q_{11} \approx 0.095$):

$$D = \sqrt{\frac{\dot{V}_n}{Q_{11} \cdot \sqrt{H_n}}} = \sqrt{\frac{24}{0.095 \cdot \sqrt{460}}} \approx 3.43\ \text{m}$$

Flow coefficient:

$$K_q = \frac{\dot{V}_n}{D^2 \cdot \sqrt{H_n}} \approx \frac{24}{3.43^2 \cdot 21.45} \approx 0.095\ \text{m}^{-1/2}\text{s}^{-1}$$

Water starting time:

$$T_w = \frac{L_p \cdot \dot{V}_n}{g \cdot A_p \cdot H_n} = \frac{500 \cdot 24}{9.81 \cdot 4\pi \cdot 460} \approx 0.44\ \text{s}$$

---

## 2. Isochronous Control Theory

An **isochronous** governor maintains the rotor speed at exactly $\omega_\text{ref}$ regardless of load by driving the steady-state speed error to zero. This is achieved by a pure integrator (no steady-state droop):

$$e(t) = \omega_\text{ref} - \omega(t)$$

$$u_\text{PID}(t) = K_p \, e(t) + K_i \int_0^t e(\tau)\,d\tau + K_d \frac{de}{dt}$$

With $K_d = 0$ and choosing gains for the Trollheim water time constant:

$$K_p = \frac{T_w}{T_w + T_{gs}} = \frac{0.44}{0.44 + 0.2} \approx 0.69 \quad \Rightarrow \text{use } K_p = 2.0$$

$$K_i = \frac{1}{T_w} \approx 2.27 \quad \Rightarrow \text{use } K_i = 2.0$$

The gate position evolves subject to rate limits (IEEE convention):

$$\dot{u}_\text{min} = -0.2\ \text{pu/s} \quad (\text{closing}), \qquad \dot{u}_\text{max} = +0.05\ \text{pu/s} \quad (\text{opening})$$

The rotor dynamics (swing equation) with turbine and generator torques:

$$J \frac{d\omega}{dt} = \tau_\text{turbine} - \tau_\text{generator} - b_v \omega$$

Mechanical power output from the Francis turbine affinity model:

$$Q = \tau_o \cdot K_q \cdot D^2 \cdot \sqrt{H}$$

$$\eta = \eta_\text{max} \left[1 - c_\eta \left(\frac{Q}{Q_\text{rated}} - 1\right)^2\right]$$

$$P_\text{mech} = \rho g Q H \eta$$

---

## 3. Julia Code

```julia
using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots

# ── Trollheim HPP — Operating Point A ────────────────────────────────────────
const H_n      = 460.0          # Nominal head [m]
const Q_n      = 24.0           # Nominal discharge [m³/s]
const P_r      = 130e6          # Rated power [W]
const n_rpm    = 375.0          # Nominal speed [RPM]
const omega_ref = 2π * n_rpm / 60  # ≈ 39.27 rad/s

# Penstock geometry (from benchmark table)
const L_p      = 500.0          # Length [m]
const D_p      = sqrt(4π / π)   # D from A_p = 4π → D = 2√π ≈ 3.545 m
const A_p      = 4π             # Cross-section [m²]

# Turbine sizing
const D_runner = 3.43           # Runner diameter [m]
const K_q      = Q_n / (D_runner^2 * sqrt(H_n))  # ≈ 0.095

# Rotor inertia — Trollheim 130 MW, M_a ≈ 8 s
const J_rotor  = 8 * P_r / omega_ref^2  # ≈ 6.7×10⁵ kg·m²

@mtkmodel TrollheimISO begin
    @components begin
        # Hydraulic circuit
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

        # Mechanical train
        rotor      = RotorInertia(
                         J       = J_rotor,
                         b_v     = 10.0,
                         omega_0 = omega_ref,
                     )
        generator  = SimpleGenerator(
                         P_rated     = P_r,
                         omega_rated = omega_ref,
                         omega_s     = omega_ref,
                         D_d         = 2.0,
                         eta_gen     = 0.98,
                     )

        # Isochronous PID governor (no droop: K_d = 0)
        governor   = PIDGovernor(
                         K_p       = 2.0,
                         K_i       = 2.0,
                         K_d       = 0.0,
                         omega_ref = omega_ref,
                         R_open    = 0.05,     # 0.05 pu/s  (benchmark)
                         R_close   = 0.20,     # 0.20 pu/s  (benchmark)
                         tau_0     = 0.8,
                     )
    end
    @equations begin
        # Hydraulic connections
        connect(reservoir.port,  penstock.port_a)
        connect(penstock.port_b, turbine.port_in)
        connect(turbine.port_out, tailwater.port)

        # Mechanical connections
        connect(turbine.shaft,          rotor.flange_turbine)
        connect(rotor.flange_generator, generator.flange)

        # Governor loop: speed → governor → gate
        connect(rotor.speed_out,   governor.speed_in)
        connect(governor.gate_out, turbine.opening)
    end
end

@named sys = TrollheimISO()
sys = structural_simplify(sys)

# Initial conditions at nominal operating point (80 % gate)
u0 = [
    sys.rotor.omega   => omega_ref,
    sys.rotor.theta   => 0.0,
    sys.governor.e_int => 0.0,
    sys.governor.tau_o => 0.8,
    sys.penstock.dm   => 1000.0 * Q_n,
]

tspan = (0.0, 30.0)
prob  = ODEProblem(sys, u0, tspan)
sol   = solve(prob, Rodas5P(), abstol=1e-6, reltol=1e-6)

# ── Plot results ─────────────────────────────────────────────────────────────
p1 = plot(sol, idxs = sys.rotor.omega,
          label = "ω [rad/s]", ylabel = "Speed (rad/s)",
          title  = "Trollheim — Isochronous Control: Rotor Speed")
hline!([omega_ref], ls=:dash, lc=:red, label="ω_ref")

p2 = plot(sol, idxs = sys.governor.tau_o,
          label = "Gate opening τₒ [-]", ylabel = "Gate opening [-]",
          title  = "Gate Position")

p3 = plot(sol, idxs = sys.turbine.P_mech ./ 1e6,
          label = "P_mech [MW]", ylabel = "Power (MW)",
          title  = "Mechanical Power Output")
hline!([P_r/1e6], ls=:dash, lc=:green, label="P_rated")

plot(p1, p2, p3, layout=(3,1), size=(800, 700), xlabel="Time (s)")
savefig("quick_examples/trollheim_isochronous.png")
println("Simulation complete. Plot saved to quick_examples/trollheim_isochronous.png")
```

---

## 4. Expected Behaviour

After a gate reduction at $t = 5\ \text{s}$:

1. **Speed rises** transiently as hydraulic power drops before the gate closes fully.
2. The **isochronous integrator** accumulates error and **drives the gate** back until $\omega = \omega_\text{ref}$ exactly.
3. **No steady-state speed deviation** — the hallmark of isochronous control.
4. **Settling time** ≈ $3\text{–}8 \times T_w \approx 1.3\text{–}3.5\ \text{s}$ depending on gains.

> **Note:** In islanded operation or single-unit grids, isochronous control is mandatory to maintain frequency. In multi-unit systems, only one unit should run isochronous — others should use droop (see Example 2).
