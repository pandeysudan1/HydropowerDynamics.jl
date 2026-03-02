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
# ---------------------------------------------------------------------------- #
@connector SignalInPort begin
    u(t), [description = "Signal (input)"]
end

@connector SignalOutPort begin
    u(t), [description = "Signal (output)"]
end
