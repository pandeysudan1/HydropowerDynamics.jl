"""
    HydroPowerDynamics.jl

Acausal, equation-based hydropower component library built on ModelingToolkit.jl.

Physical domains covered
  • Hydraulic  (pressure / mass-flow)
  • Mechanical rotational (angular velocity / torque)
  • Signal  (causal scalar)

Component families
  • Hydraulic  : Reservoir, Penstock, SurgeTank, DraftTube, GuideVane
  • Turbine    : FrancisTurbineAffinity, PeltonTurbine
  • Mechanical : RotorInertia, SimpleGenerator
  • Governor   : PIDGovernor, GGOV1Governor

Quick-start
-----------
```julia
using HydroPowerDynamics, ModelingToolkit, OrdinaryDiffEq
# See test/ for worked examples
```
"""
module HydroPowerDynamics

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import OrdinaryDiffEq
using DataInterpolations

# Re-export MTK fundamentals so users need only `using HydroPowerDynamics`
export t, D

# ── sub-modules / source files ───────────────────────────────────────────────
include("connectors.jl")
include("hydraulic.jl")
include("turbine.jl")
include("mechanical.jl")
include("governor.jl")
include("utils.jl")

# ── connectors ───────────────────────────────────────────────────────────────
export HydraulicPort, RotationalPort, SignalInPort, SignalOutPort

# ── hydraulic components ─────────────────────────────────────────────────────
export Reservoir, Penstock, SurgeTank, DraftTube, GuideVane

# ── turbine models ───────────────────────────────────────────────────────────
export FrancisTurbineAffinity, PeltonTurbine

# ── mechanical / generator ───────────────────────────────────────────────────
export RotorInertia, SimpleGenerator

# ── governor / control ───────────────────────────────────────────────────────
export PIDGovernor, GGOV1Governor

# ── system utilities ─────────────────────────────────────────────────────────
export gross_head, net_head, darcy_head_loss, power_cascade, joukowsky_pressure, wave_speed
export critical_closure_time, unit_speed, unit_discharge, plant_efficiency, hydraulic_efficiency

end # module HydroPowerDynamics
