# Surge Tank — Reference Document

## 1. Real-Life Introduction

### The water hammer problem

When a turbine governor closes the guide vane rapidly (e.g. during load rejection
or emergency shutdown), the large mass of water in the headrace tunnel cannot stop
instantly. The kinetic energy of the water column must go somewhere. Without a
pressure relief mechanism:

1. The column decelerates and a compression wave travels upstream at the **acoustic
   wave speed** $a \approx 1200\ \text{m/s}$, generating an overpressure:

   $$\Delta p_\text{Joukowsky} = \rho \, a \, \Delta v$$

   For Trollheim ($\rho=1000$, $a \approx 1200\ \text{m/s}$, $\Delta v \approx Q/A \approx 37/31 \approx 1.2\ \text{m/s}$):

   $$\Delta p \approx 1000 \times 1200 \times 1.2 = 1.44\ \text{MPa} \approx 147\ \text{m head}$$

   This is **40 % of the 371 m nominal head** — large enough to rupture a penstock.

2. If the closure is **slow** (longer than the critical reflection time $T_c = 2L/a$),
   high-frequency water hammer is avoided. But a **low-frequency surge** still builds
   up in the headrace over a period of 60–120 s, oscillating back and forth and
   causing the turbine head to fluctuate.

### Why a surge tank helps

A surge tank is an open-top (or pressurised) chamber placed at the junction of the
headrace tunnel and the penstock. It solves both problems simultaneously:

| Without surge tank | With surge tank |
|---|---|
| Full Joukowsky pressure acts on penstock | Tank absorbs flow excess → penstock sees only quasi-static ΔH |
| Tunnel inertia interacts directly with turbine | Tunnel inertia decoupled by free surface |
| Governor must fight head oscillations | Head at turbine ≈ tank level (slowly varying) |

In **grid-connected** operation the surge tank is critical for **governor stability**:
the hydraulic system has a negative initial response (increasing gate → flow lags →
temporary head drop → power drops before rising). The surge tank reduces the water
starting time $T_w$ seen by the governor, widening the stable gain range.

### Trollheim layout

```
  Reservoir (H_res = 384.9 m mean SS)
      │
      │  Headrace tunnel: L_t = 4 496 m, D_t = 6.3 m
      │  (intake race #1 81.5 m + #2 395 m + #3 4020 m from benchmark table)
      │
    [Surge tank]  D = 3.4 m → A_t ≈ 9.1 m²,  height extent 75.5 m
      │                        riser shaft: D_r, L_r, ε_r  (tunable)
      │
      │  Penstock: L_p = 508 m, D_p ≈ 3.9 m (avg of two sections)
      │
    [Francis turbine + generator]
      │
  Tailwater (H_tail ≈ 0 m datum)
```

---

## 2. Mathematical Construction

### 2.1 State variables

| Symbol | Description | Units |
|---|---|---|
| $Z(t)$ | Water surface elevation in surge tank | m |
| $Q_t(t)$ | Volumetric flow in headrace tunnel | m³/s |
| $Q(t)$ | Turbine flow (algebraic from FA.1) | m³/s |
| $\omega(t)$ | Rotor angular velocity | rad/s |

### 2.2 Headrace tunnel momentum equation

Applying Newton's second law to the tunnel water column (length $L_t$,
cross-section $A_t$, ignoring compressibility):

$$\frac{L_t}{g A_t} \frac{dQ_t}{dt} = H_\text{res} - Z(t) - h_{f,t}(Q_t) \tag{TU.1}$$

where the friction head loss in the tunnel is:

$$h_{f,t} = f_D(Re_t,\, D_t,\, \varepsilon_t) \cdot \frac{L_t}{D_t} \cdot \frac{Q_t^2}{2 g A_t^2} \tag{TU.2}$$

This is the **Darcy-Weisbach** equation with a dynamically computed friction factor
--- see Section 3 for the multi-regime implementation.

### 2.3 Surge tank continuity

The free surface rises when tunnel inflow exceeds turbine demand:

$$A_s \frac{dZ}{dt} = Q_t - Q_\text{turb} \tag{ST.1}$$

where $A_s$ is the **tank cross-sectional area** (not the riser shaft area).

### 2.4 Riser throttle friction

A real surge tank connects to the tunnel through a **riser shaft** (sometimes a
separate throttle orifice) of diameter $D_r$, length $L_r$, and roughness $\varepsilon_r$.
The riser provides a flow-dependent pressure drop that:
- **damps** surge oscillations (adds energy dissipation)
- Is the **primary tuning parameter** for surge tank behaviour

The riser friction head loss (Darcy-Weisbach, sign-preserving):

$$h_{f,r} = f_D(Re_r,\, D_r,\, \varepsilon_r) \cdot \frac{L_r}{D_r} \cdot \frac{Q_r\,|Q_r|}{2 g A_r^2} \tag{ST.3}$$

The pressure at the hydraulic port thus becomes:

$$p_\text{port} = p_\text{atm} + \rho g Z - \text{sign}(Q_r) \cdot \rho g\, h_{f,r} \tag{ST.2}$$

A higher roughness $\varepsilon_r$ → larger $f_D$ → stronger damping of surge oscillations,
but also greater steady-state head loss.

### 2.5 Effective turbine head

With surge tank present, the turbine sees:

$$H_\text{eff}(t) = Z(t) - h_{f,p}(Q_t) - H_\text{tail} \tag{TU.3}$$

where $h_{f,p}$ is the penstock friction head loss. Without a surge tank,
$H_\text{eff} = H_\text{res} - h_{f,t} - h_{f,p}$ (constant at rated flow).

### 2.6 Surge tank natural period

The un-damped natural period of the headrace-tank oscillation is:

$$T_s = 2\pi \sqrt{\frac{L_t \, A_s}{g \, A_t^2}} \tag{ST.4}$$

For Trollheim ($L_t = 4496\ \text{m}$, $A_s = 9.08\ \text{m}^2$, $A_t = 31.2\ \text{m}^2$):

$$T_s = 2\pi \sqrt{\frac{4496 \times 9.08}{9.81 \times 31.2^2}}
      = 2\pi \sqrt{\frac{40\,823}{9\,541}} \approx 2\pi \times 2.07 \approx 72.6\ \text{s}$$

This slow oscillation (≈ 73 s period) interacts with the governor response
(settling time ≈ 20–40 s) and can cause sustained power swings if the governor
is not de-tuned or if riser damping is insufficient.

### 2.7 Stability criterion (Thoma criterion)

The minimum surge tank area for stable regulation is (Thoma 1910):

$$A_s \geq A_{s,\text{Thoma}} = \frac{L_t \, A_t}{2 \, h_{f,t,0} \, T_m} \tag{ST.5}$$

where $h_{f,t,0}$ is the rated friction head loss in the tunnel and $T_m = J \omega^2 / P_r$
is the mechanical starting time. Tanks smaller than $A_{s,\text{Thoma}}$ lead to
growing oscillations under isochronous control.

---

## 3. Friction Factor: Multi-Regime Darcy-Weisbach

### 3.1 Why a single turbulent formula is insufficient

The Swamee-Jain approximation (used in the original `Penstock` model):

$$f_D^\text{SJ} = \frac{0.25}{\left[\log_{10}\!\left(\frac{\varepsilon}{3.7D} + \frac{5.74}{Re^{0.9}}\right)\right]^2}$$

diverges for $Re \to 0$ (gives $f_D \to \infty$). During startup, shutdown, and
load-rejection transients, flow passes through the **laminar regime** ($Re < 2100$)
where the exact Hagen-Poiseuille result applies: $f = 64/Re$.

### 3.2 Three-regime implementation

The `darcy_factor(Re, D, ε)` function in `utils.jl` implements:

$$f_D(Re) = \begin{cases}
\dfrac{64}{Re} & Re \leq 2100 \quad \text{(laminar, Hagen-Poiseuille)} \\[8pt]
\text{cubic Hermite spline} & 2100 < Re < 2300 \quad \text{(smooth transition)} \\[8pt]
\dfrac{1}{\left[2\log_{10}\!\left(\dfrac{\varepsilon}{3.7D} + \dfrac{5.7}{Re^{0.9}}\right)\right]^2} & Re \geq 2300 \quad \text{(turbulent, Swamee-Jain variant)}
\end{cases}$$

The cubic Hermite spline in the transition zone matches **value and first derivative**
at both boundaries ($Re = 2100$ and $Re = 2300$), giving a $C^1$-continuous function.
This is important for ODE solver Jacobian accuracy — a discontinuous $df_D/dRe$ at
the laminar-turbulent boundary can cause the solver to take unnecessarily small steps
or fail near $Re \approx 2100$.

### 3.3 Effect of roughness height on surge tank damping

The riser roughness $\varepsilon_r$ controls the level of turbulent friction in the
tank connection. Physically meaningful values for different construction types:

| Lining | $\varepsilon_r$ [m] | Expected effect |
|---|---|---|
| Smooth steel pipe | $4.6 \times 10^{-5}$ | Minimal damping, slow decay |
| Concrete / shotcrete | $1.0 \times 10^{-3}$ | Moderate damping, ~3–5 cycles |
| Unlined blasted rock | $1.0 \times 10^{-2}$ | Heavy damping, 1–2 cycles |

Trollheim's surge tank shaft (D = 3.4 m, excavated rock) likely falls in the
**unlined rock** to **shotcrete** range. The `sim_test` script sweeps all three.

---

## 4. Code Implementation in `HydroPowerDynamics.jl`

### 4.1 `darcy_factor` in `src/utils.jl`

```julia
@register_symbolic darcy_factor(Re, D, e)       # expose to MTK symbolic engine

function darcy_factor(Re, D, e)
    Re_l, Re_t = 2100.0, 2300.0
    f_l(r) = 64.0 / r
    f_t(r) = 1.0 / (2.0 * log10(e / (3.7D) + 5.7 / (r^0.9)))^2

    Re <= 0    && return 0.0
    Re <= Re_l && return f_l(Re)
    Re >= Re_t && return f_t(Re)

    # Cubic Hermite spline in [Re_l, Re_t]
    y1, y2  = f_l(Re_l), f_t(Re_t)
    dy1     = -64.0 / Re_l^2
    δ       = 1e-3 * Re_t
    dy2     = (f_t(Re_t + δ) - f_t(Re_t - δ)) / (2δ)
    X       = [Re_l^3 Re_l^2 Re_l 1.0
               Re_t^3 Re_t^2 Re_t 1.0
               3Re_l^2 2Re_l 1.0 0.0
               3Re_t^2 2Re_t 1.0 0.0]
    k       = X \ [y1, y2, dy1, dy2]
    return k[1]*Re^3 + k[2]*Re^2 + k[3]*Re + k[4]
end
```

The `@register_symbolic` macro tells ModelingToolkit to treat `darcy_factor` as an
opaque callable during symbolic simplification and use automatic differentiation
(ForwardDiff) for the Jacobian — no symbolic expansion of the conditional logic.

### 4.2 `Penstock` update in `src/hydraulic.jl`

The single-regime Swamee-Jain line is replaced by `darcy_factor`:

```julia
# OLD — wrong for Re < 2100
f_D ~ 0.25 / (log10(e / (3.7 * D) + 5.74 / (Re^0.9 + 1e-10)))^2

# NEW — correct across all regimes
f_D ~ darcy_factor(Re, D, e)
```

### 4.3 `SurgeTank` update in `src/hydraulic.jl`

Three new parameters control the riser shaft:

```julia
D_riser = 3.4,     # riser diameter [m]   (Trollheim default from geometry table)
L_riser = 87.0,    # riser length [m]     (height extent from benchmark)
e_riser = 1.0e-2,  # roughness [m]        (unlined rock default)
```

New variables: `Re_riser(t)`, `f_riser(t)`, `v_riser(t)`, `dp_riser(t)`.

The port pressure becomes:

```julia
port.p ~ p_atm + rho * g * Z
         - sign(port.dm + 1e-20) * dp_riser   # sign-preserving riser friction
```

### 4.4 Flat model in `sim_test/iso_surgetank_comparison.jl`

Because `HydroPowerDynamics.jl` follows the flat `@mtkmodel` pattern (no
`connect()` wiring), the surge tank scenario adds two new state variables to the
existing ISO model:

| New state | Equation | Source |
|---|---|---|
| `Z(t)` — tank level | `A_s * D(Z) ~ Q_t - Q_turb` | ST.1 |
| `Q_t(t)` — tunnel flow | `(L_t/(g*A_tun)) * D(Q_t) ~ H_res - Z - h_ft` | TU.1 |

Effective turbine head: `H_eff = Z + Z_ref - h_fp - H_tail`
where `Z_ref` converts tank level to absolute head above the turbine datum.

---

## 5. References

1. **Wylie & Streeter** (1993), *Fluid Transients in Systems*, Prentice Hall.
2. **Thoma, D.** (1910), Zur Theorie des Wasserschlosses bei selbsttätig geregelten Turbinenanlagen, *Oldenbourg, Munich*.
3. **IEC 60193:2019** — Hydraulic turbines, storage pumps and pump-turbines — Model acceptance tests.
4. **IEEE Std 1207-2011** — Guide for the Application of Turbine Governing Systems.
5. **Kundur, P.** (1994), *Power System Stability and Control*, Ch. 15 (Hydraulic turbine governing).
6. Trollheim HPP waterway geometry: Energies 2019 figure/table (compiled in `benchmarks/trollheim_hydropower.md`).
