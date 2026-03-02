# test_utils.jl – Unit tests for HydroPowerDynamics utility functions (no MTK required)
using Test
using HydroPowerDynamics

@testset "Utility Functions" begin

    @testset "Head equations (Section 6.1)" begin
        H_g = gross_head(150.0, 5.0)
        @test H_g ≈ 145.0

        H_n = net_head(145.0, 3.5, 1.0)
        @test H_n ≈ 140.5

        # Darcy head loss: f=0.02, L=500, D=2, v=2 m/s, g=9.81
        hf = darcy_head_loss(0.02, 500.0, 2.0, 2.0)
        @test hf ≈ 0.02 * (500/2) * (4 / (2*9.81))  rtol = 1e-6
    end

    @testset "Power cascade (Section 6.2)" begin
        pc = power_cascade(1000.0, 9.81, 80.0, 150.0, 145.0, 0.91, 0.97)
        @test pc.P_hydraulic  ≈ 1000 * 9.81 * 80 * 150
        @test pc.P_available  ≈ 1000 * 9.81 * 80 * 145
        @test pc.P_turbine    ≈ pc.P_available * 0.91  rtol = 1e-9
        @test pc.P_electrical ≈ pc.P_turbine   * 0.97  rtol = 1e-9
        @test pc.eta_plant    < 1.0
        @test pc.eta_plant    > 0.0
    end

    @testset "Plant efficiency (Eq. PC.7)" begin
        @test plant_efficiency(0.91, 0.99, 0.97, 0.995) ≈ 0.91 * 0.99 * 0.97 * 0.995
    end

    @testset "Water hammer (Section 6.3)" begin
        # Steel pipe: E = 200 GPa, e_wall = 0.02 m
        a = wave_speed(2.0e9, 1000.0, 2.0, 200e9, 0.02)
        @test 900.0 < a < 1500.0   # typical range for steel penstock

        dp = joukowsky_pressure(1000.0, a, 2.0)
        @test dp > 0.0
        @test dp ≈ 1000.0 * a * 2.0

        Tc = critical_closure_time(500.0, a)
        @test Tc ≈ 2 * 500.0 / a
    end

    @testset "Dimensionless turbine parameters (Section 3.1)" begin
        n11 = unit_speed(300.0, 2.5, 100.0)
        @test n11 ≈ 300.0 * 2.5 / sqrt(100.0)   # = 75 rpm·m^0.5

        Q11 = unit_discharge(80.0, 2.5, 100.0)
        @test Q11 ≈ 80.0 / (2.5^2 * sqrt(100.0))

        eta = hydraulic_efficiency(40e6, 1000.0, 9.81, 80.0, 100.0)
        @test 0.0 < eta < 1.0
        @test eta ≈ 40e6 / (1000 * 9.81 * 80 * 100)
    end
end
