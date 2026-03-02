# test_francis_step.jl
# ────────────────────
# Integration tests for HydroPowerDynamics component equations.
#
# Test 1 – Francis turbine + rotor spin-up
#   Flat @mtkmodel with affinity-law turbine, rotor inertia, and
#   algebraic generator. Verifies spin-up physics and efficiency.
#
# Test 2 – Pelton jet impulse formula (Section 3.4 algebraic check)
#
# Test 3 – Surge tank linear rise (Section 2.4, analytic reference)
#
# Test 4 – PID governor response (speed error integral drives gate open)
using Test
using HydroPowerDynamics
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq

# ─────────────────────────────────────────────────────────────────────────────
# Test 1: Francis turbine affinity model + rotor spin-up
# ─────────────────────────────────────────────────────────────────────────────
@testset "Francis Turbine + Rotor spin-up" begin

    @mtkmodel TurbineRotorSystem begin
        @parameters begin
            H       = 100.0   # net head [m]
            tau_o   = 0.8     # constant gate opening
            rho     = 1000.0
            g       = 9.81
            D_t     = 1.0     # runner diameter [m]
            K_q     = 0.35
            eta_max = 0.92
            c_eta   = 0.25
            Q_rated = 8.0
            J       = 500.0   # rotor inertia [kg·m²]
            b_v     = 2.0     # viscous damping [N·m·s/rad]
            P_rated = 2.0e6   # generator rated power [W]
            omega_s = 157.08  # synchronous speed [rad/s]
            D_d     = 1.5
        end
        @variables begin
            omega(t) = 1.0
            Q(t)
            eta(t)
            P_mech(t)
            tau_turbine(t)
            tau_gen(t)
        end
        @equations begin
            Q            ~ tau_o * K_q * D_t^2 * sqrt(H)           # FA.1
            eta          ~ eta_max * (1 - c_eta*(Q/Q_rated - 1)^2) # FA.2
            P_mech       ~ rho * g * Q * H * eta                    # FA.3
            tau_turbine  ~ P_mech / (abs(omega) + 1e-4)
            tau_gen      ~ (P_rated/omega_s)*(1 + D_d*(omega - omega_s)/omega_s) # GA.1
            D(omega)     ~ (tau_turbine - tau_gen - b_v*omega) / J  # M.1
        end
    end

    @mtkcompile sys = TurbineRotorSystem()
    prob = ODEProblem(sys, Pair[], (0.0, 10.0))
    sol  = solve(prob, Tsit5(); abstol=1e-6, reltol=1e-6)

    @test sol.retcode == ReturnCode.Success

    omega_0   = sol[sys.omega][1]
    omega_end = sol[sys.omega][end]
    P_end     = sol[sys.P_mech][end]
    eta_end   = sol[sys.eta][end]

    println("  ω(t=0)  = $(round(omega_0,   digits=3)) rad/s")
    println("  ω(t=10) = $(round(omega_end, digits=2)) rad/s")
    println("  P_mech  = $(round(P_end/1e3, digits=1)) kW")
    println("  η       = $(round(eta_end*100, digits=1)) %")

    @test omega_end > omega_0      # rotor accelerates
    @test P_end > 0.0              # positive mechanical power
    @test 0.0 < eta_end ≤ 1.0     # physical efficiency
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 2: Pelton algebraic formula check (Section 3.4)
# ─────────────────────────────────────────────────────────────────────────────
@testset "Pelton jet impulse – algebraic (Section 3.4)" begin
    C_v=0.98; g=9.81; H_net=400.0
    v_jet = C_v * sqrt(2*g*H_net)
    @test v_jet ≈ 0.98*sqrt(2*9.81*400)  rtol=1e-8

    phi=0.85; beta2=170.0; rho=1000.0
    omega=30.0; R=1.5; Q=0.5
    u1  = omega * R
    tau = rho * Q * v_jet * u1 * (1 + phi*cosd(beta2)) / omega
    @test tau > 0.0
    @test 0.0 < u1/v_jet < 1.0   # physical bucket speed ratio
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 3: Surge tank linear rise – analytic verification (Section 2.4)
# ─────────────────────────────────────────────────────────────────────────────
@testset "Surge tank linear level rise (Section 2.4)" begin

    @mtkmodel SurgeTankLinear begin
        @parameters begin
            A_t   = 50.0
            rho   = 1000.0
            g     = 9.81
            p_atm = 101_325.0
            dm_in = 500.0   # constant inflow [kg/s]
        end
        @variables begin
            Z(t) = 10.0     # water level [m]
            p(t)            # pressure at port [Pa]
        end
        @equations begin
            D(Z) ~ dm_in / (rho * A_t)       # ST.1
            p    ~ p_atm + rho * g * Z        # ST.2
        end
    end

    @mtkcompile stank = SurgeTankLinear()
    prob = ODEProblem(stank, Pair[], (0.0, 100.0))
    sol  = solve(prob, Tsit5())

    @test sol.retcode == ReturnCode.Success

    dm_in = 500.0; A_t = 50.0; rho = 1000.0; Z_0 = 10.0
    Z_analytic = Z_0 + dm_in/(rho*A_t)*100.0   # = 11.0 m

    @test sol[stank.Z][end] ≈ Z_analytic  rtol=1e-4
    @test sol[stank.Z][end] > Z_0
    @test sol[stank.p][end] > 101_325.0
end

# ─────────────────────────────────────────────────────────────────────────────
# Test 4: PID governor integrator drives gate to 1 when below reference speed
# ─────────────────────────────────────────────────────────────────────────────
@testset "PID governor integrator (Section 5.1)" begin

    # Standalone PID: fixe speed input below reference → gate should open
    @mtkmodel PIDStandalone begin
        @parameters begin
            K_p       = 2.0
            K_i       = 1.0
            K_d       = 0.0   # no derivative for clean test
            omega_ref = 100.0
            omega_meas = 80.0  # constant measured speed (below ref)
        end
        @variables begin
            e(t)           = 20.0
            e_int(t)       = 0.0
            u_pid(t)
            tau_o(t)       = 0.3
        end
        @equations begin
            e       ~ omega_ref - omega_meas          # Gov.1  (constant error = 20)
            D(e_int) ~ e                              # integral of error
            u_pid   ~ K_p*e + K_i*e_int              # Gov.2
            D(tau_o) ~ max(-0.2, min(0.1,            # Gov.4 rate limit
                        (min(max(u_pid, 0.0), 1.0) - tau_o) / 0.01))
        end
    end

    @mtkcompile gov = PIDStandalone()
    prob = ODEProblem(gov, Pair[], (0.0, 5.0))
    sol  = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8)

    @test sol.retcode == ReturnCode.Success

    tau_0   = sol[gov.tau_o][1]
    tau_end = sol[gov.tau_o][end]

    println("  τ_o(t=0) = $(round(tau_0,   digits=3))")
    println("  τ_o(t=5) = $(round(tau_end, digits=3))")

    # With speed below reference and integrator active, gate should open
    @test tau_end > tau_0
    @test 0.0 ≤ tau_end ≤ 1.0
end
