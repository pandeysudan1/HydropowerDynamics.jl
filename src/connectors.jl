# ============================================================================ #
# connectors.jl  –  Port definitions for all physical domains
#
# Each connector defines:
#   • across variable  – equalised at every connect() node    (potential)
#   • through variable – sums to zero at every connect() node (flow)
#
# The `connect = Flow` attribute tells MTK to generate the Kirchhoff sum
# equations automatically.
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Hydraulic domain
# across : p  [Pa]        pressure
# through: dm [kg/s]      mass flow rate  (positive → entering component)
# ---------------------------------------------------------------------------- #
@connector HydraulicPort begin
    p(t),  [description = "Pressure [Pa]"]
    dm(t), [description = "Mass flow rate [kg/s]", connect = Flow]
end

# ---------------------------------------------------------------------------- #
# Mechanical rotational domain
# across : phi [rad], omega [rad/s]
# through: tau [N·m]   torque          (positive → entering component)
# ---------------------------------------------------------------------------- #
@connector RotationalPort begin
    phi(t),   [description = "Angular position [rad]"]
    omega(t), [description = "Angular velocity [rad/s]"]
    tau(t),   [description = "Torque [N·m]",  connect = Flow]
end

# ---------------------------------------------------------------------------- #
# Signal domain  (causal – no flow variable, no Kirchhoff sum)
#
# Delegated to ModelingToolkitStandardLibrary.Blocks so that:
#   • connect(sensor.w, :y, governor.speed_in)  — analysis-point tagging works
#   • Blocks.get_sensitivity / get_looptransfer  — can intercept these nodes
#   • All MTSL causal blocks (Feedback, FirstOrder, LimPID…) are plug-compatible
#
# Blocks.RealInput  has port variable  u(t)  — same name as our old custom port,
# so all existing components (turbine, mechanical) are drop-in compatible.
# ---------------------------------------------------------------------------- #
const SignalInPort  = Blocks.RealInput    # causal input  connector  (variable: u)
const SignalOutPort = Blocks.RealOutput   # causal output connector  (variable: u)
