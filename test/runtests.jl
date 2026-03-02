# test/runtests.jl – HydroPowerDynamics test suite
using Test

@testset "HydroPowerDynamics" begin
    include("test_utils.jl")
    include("test_francis_step.jl")
end
