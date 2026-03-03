# Plan: MTSL-Style Architecture with Analysis Points

**Goal**: Rebuild `connection_test.jl` using proper sensor components and MTK analysis
points so that `ControlSystemsBase` can compute `S(s)`, `T(s)`, `L(s)` directly from
the acausal hydropower model.

---

## Current Problem (connection_test.jl)

```
Speed feedback is an algebraic equation hack:
    governor.speed_in.u ~ rotor.omega          ← no analysis point possible
    governor.power_in.u ~ turbine.P_mech       ← same issue

Gate signal:
    turbine.opening.u   ~ governor.gate_out.u  ← no analysis point possible
```

Because there is no intermediate connector node, there is nowhere to insert `:y`
or `:u` tags, so `Blocks.get_sensitivity(model, :y)` cannot be called.

---

## Target Architecture (DC Motor pattern applied to hydropower)

```
                 ┌─────────────────────────────────────────────┐
                 │         TrollheimGGOV1 (@mtkmodel)           │
                 │                                              │
  reservoir ─── penstock ─── turbine ─── speed_sensor ──:y──▶ governor
                                │                               │
                              shaft ──── rotor ──── generator   │
                                │                               │
                         power_sensor ──────────────────────▶  │
                                                                │
                                          :u ◀── gate_out ─────┘
                                           │
                                        turbine.opening
                 └─────────────────────────────────────────────┘
```

Key loops:
- **Analysis point `:y_omega`** — between `speed_sensor.w` and `governor.speed_in`
- **Analysis point `:u_gate`**  — between `governor.gate_out` and `turbine.opening`

---

## Work Packages

### WP-1  Add sensor components to `src/mechanical.jl`

#### 1a. `RotationalSpeedSensor`
| | |
|---|---|
| Ports | `flange::RotationalPort`, `w::SignalOutPort` |
| Equations | `flange.tau ~ 0` (ideal — draws no energy), `w.u ~ flange.omega` |
| Connection | `connect(rotor.flange_generator, speed_sensor.flange)` |
| Maps to | MTSL `SpeedSensor` |

#### 1b. `MechanicalPowerSensor`
| | |
|---|---|
| Ports | `flange_a::RotationalPort`, `flange_b::RotationalPort`, `P::SignalOutPort` |
| Equations | Pass-through: `flange_b.tau ~ -flange_a.tau`, `flange_b.omega ~ flange_a.omega` |
| | Output: `P.u ~ flange_a.tau * flange_a.omega` |
| Connection | Inserted in-line: `turbine.shaft ── power_sensor ── rotor.flange_turbine` |
| Maps to | MTSL `PowerSensor` (in-line, zero power loss) |

Exports to add in `HydroPowerDynamics.jl`:
```julia
export RotationalSpeedSensor, MechanicalPowerSensor
```

---

### WP-2  Verify signal connector compatibility

`SignalInPort` / `SignalOutPort` both declare only `u(t)` — same as MTSL's
`RealInput` / `RealOutput`. MTK analysis points work on any such causal connector
pair, so **no changes to `connectors.jl` are needed**.

Confirm the tag syntax works:
```julia
connect(speed_sensor.w,        :y_omega,  governor.speed_in)
connect(governor.gate_out,     :u_gate,   turbine.opening)
```
These generate `governor.speed_in.u ~ speed_sensor.w.u` with an interceptable
intermediate node named `:y_omega` / `:u_gate`.

---

### WP-3  New simulation file: `connection_test_v2.jl`

Location: `quick_examples/sim_test/connection_test_v2.jl`

#### 3a. Component wiring changes vs. current `connection_test.jl`

| Signal | Old (algebraic hack) | New (sensor + connect) |
|---|---|---|
| Speed | `governor.speed_in.u ~ rotor.omega` | `connect(rotor.flange_gen, speed_sensor.flange)` + `connect(speed_sensor.w, :y_omega, governor.speed_in)` |
| Power | `governor.power_in.u ~ turbine.P_mech` | `connect(turbine.shaft, power_sensor.flange_a)` + `connect(power_sensor.flange_b, rotor.flange_turbine)` + `connect(power_sensor.P, governor.power_in)` |
| Gate | `turbine.opening.u ~ governor.gate_out.u` | `connect(governor.gate_out, :u_gate, turbine.opening)` |

#### 3b. Additional imports needed
```julia
using ModelingToolkitStandardLibrary.Blocks   # for get_sensitivity / get_looptransfer
using ControlSystemsBase                       # for ss(), bodeplot(), nyquistplot()
```

#### 3c. Full `@equations` block (new signals section)
```julia
@equations begin
    # ── Physical domains (unchanged) ────────────────────────────
    connect(reservoir.port,   penstock.port_a)
    connect(penstock.port_b,  power_sensor.flange_a)   # power sensor in-line
    connect(power_sensor.flange_b, turbine.port_in)    # NOTE: mechanical, not hydraulic
    # Actually power sensor goes on mechanical shaft, not hydraulic:
    connect(penstock.port_b,  turbine.port_in)          # hydraulic unchanged
    connect(turbine.port_out, tailwater.port)
    connect(turbine.shaft,    power_sensor.flange_a)    # mechanical power tap
    connect(power_sensor.flange_b, rotor.flange_turbine)
    connect(rotor.flange_generator, speed_sensor.flange)
    connect(rotor.flange_generator, generator.flange)

    # ── Cross-domain bridges with analysis points ────────────────
    connect(speed_sensor.w,    :y_omega,  governor.speed_in)   # plant output tag
    connect(governor.gate_out, :u_gate,   turbine.opening)     # plant input tag
    connect(power_sensor.P,               governor.power_in)   # no tag needed (inner signal)
end
```

#### 3d. Control analysis block (after simulation)
```julia
# ── Linearise around steady-state operating point ────────────────
op = Dict(unknowns(sys) .=> 0.0)    # or use solved steady-state values

# Sensitivity function  S(s) = 1 / (1 + P·C)
matrices_S, _ = Blocks.get_sensitivity(model, :y_omega; op)
S = ss(matrices_S...) |> minreal

# Complementary sensitivity  T(s) = P·C / (1 + P·C)
matrices_T, _ = Blocks.get_comp_sensitivity(model, :y_omega; op)
T_cs = ss(matrices_T...)

# Loop transfer function  L(s) = P·C  (negative feedback sign corrected)
matrices_L, _ = Blocks.get_looptransfer(model, :y_omega; op)
L = -ss(matrices_L...)

# Gain margin, phase margin
gm, pm, wgc, wpc = margin(L)

# Plots
p_bode   = bodeplot([S, T_cs]; label=["S(s)" "T(s)"], plotphase=false,
                    title="Trollheim GGOV1 — Sensitivity functions")
p_nyq    = nyquistplot(L; label="L(s)",
                       title="Trollheim GGOV1 — Nyquist")
```

---

### WP-4  Update `Project.toml` compat

`ControlSystemsBase` is already in deps. Ensure version compat is set:
```toml
[compat]
ControlSystemsBase = "1.1"
```
Already present. ✅

`ModelingToolkitStandardLibrary = "2"` — add `using` in the script only (not a
module dependency of `HydroPowerDynamics.jl` itself; it's a simulation-level
tool).

---

## Implementation Order

| Step | File(s) | Action |
|------|---------|--------|
| 1 | `src/mechanical.jl` | Add `RotationalSpeedSensor` |
| 2 | `src/mechanical.jl` | Add `MechanicalPowerSensor` |
| 3 | `src/HydroPowerDynamics.jl` | Export new components |
| 4 | `quick_examples/sim_test/connection_test_v2.jl` | New script with sensor wiring + `:y_omega` / `:u_gate` |
| 5 | same file | Add linearisation + Bode/Nyquist block |
| 6 | same file | Run and verify `sol.retcode == Success` |
| 7 | same file | Run sensitivity analysis, save plots |

Total new/changed files: **3**  
New components: **2** (`RotationalSpeedSensor`, `MechanicalPowerSensor`)

---

## Expected Outputs

```
quick_examples/sim_test/images/
    connection_test_v2.png          ← time-domain: ω, gate, P_mech, Δω
    dc_frequency_response.png       ← Bode: S(s) and T(s) with GGOV1 droop
    nyquist_looptransfer.png        ← Nyquist: L(s), stability margins printed
```

Terminal output will include:
```
Gain margin  GM = X.X dB   at ω = X.XX rad/s
Phase margin PM = XX.X°    at ω = X.XX rad/s
```

---

## Why this matters for Trollheim

The ISO governor and GGOV1 both have tunable gains. With `GM > 6 dB` and
`PM > 45°` as stability targets, the Bode plot immediately shows whether the
current `sigma=0.10`, `T_gov=1.75 s`, `K_i=5.0` combination is robust — without
manually perturbing the nonlinear ODE and fitting frequency responses by hand.
