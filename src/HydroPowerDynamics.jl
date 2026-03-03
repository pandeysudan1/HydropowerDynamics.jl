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

# ModelingToolkitStandardLibrary — used as base for all causal signal blocks.
# Blocks.RealInput / RealOutput become the canonical signal connector types.
# All governor and sensor components are built from Blocks primitives.
using ModelingToolkitStandardLibrary.Blocks

# Re-export MTK fundamentals + Blocks so users need only `using HydroPowerDynamics`
export t, D
export Blocks   # expose full MTSL Blocks namespace: Feedback, FirstOrder, LimPID…

# ── sub-modules / source files ───────────────────────────────────────────────
include("utils.jl")       # must be first: defines darcy_factor used by hydraulic.jl
include("connectors.jl")
include("hydraulic.jl")
include("turbine.jl")
include("mechanical.jl")
include("governor.jl")

# ── connectors ───────────────────────────────────────────────────────────────
# SignalInPort / SignalOutPort are now aliases for Blocks.RealInput / RealOutput.
# Existing component code (turbine, mechanical) is unchanged — the aliases
# provide full backward compatibility while enabling analysis-point tagging.
export HydraulicPort, RotationalPort, SignalInPort, SignalOutPort

# ── hydraulic components ─────────────────────────────────────────────────────
export Reservoir, Penstock, SurgeTank, DraftTube, GuideVane

# ── turbine models ───────────────────────────────────────────────────────────
export FrancisTurbineAffinity, PeltonTurbine

# ── mechanical / generator ───────────────────────────────────────────────────
export RotorInertia, SimpleGenerator, RotationalSpeedSensor, MechanicalPowerSensor

# ── governor / control ───────────────────────────────────────────────────────
export PIDGovernor, GGOV1Governor

# ── system utilities ─────────────────────────────────────────────────────────
export gross_head, net_head, darcy_head_loss, power_cascade, joukowsky_pressure, wave_speed
export critical_closure_time, unit_speed, unit_discharge, plant_efficiency, hydraulic_efficiency
export darcy_factor, swamee_jain

end # module HydroPowerDynamics
