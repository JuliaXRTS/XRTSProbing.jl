module ThomsonTests

using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

include("groundtruths/thomson_scattering.jl")


const OMEGAS = (rand(RNG), 1e2 * rand(RNG), rand(RNG), 1e3 * rand(RNG), 1e4 * rand(RNG))
const CTHS = (-1.0, 0.0, 1.0, rand(RNG), -rand(RNG))
const PHIS = [0, 2 * pi, rand(RNG) * 2 * pi]

const PROC = ThomsonScattering()
const MODEL = PerturbativeQED()
const PSDEF = ElabPhotonSphSystem()

@testset "static properties" begin

    @test incoming_particles(PROC) == (Electron(), Photon())
    @test outgoing_particles(PROC) == (Electron(), Photon())

    @test QEDbase.in_phase_space_dimension(PROC, MODEL) == 1
    @test QEDbase.out_phase_space_dimension(PROC, MODEL) == 2

end

@testset "($om,$cth,$phi)" for (om, cth, phi) in Iterators.product(OMEGAS, CTHS, PHIS)
    in_ps = (om,)
    out_ps = (cth, phi)
    psp = PhaseSpacePoint(PROC, MODEL, PSDEF, in_ps, out_ps)

    @testset "particle momenta" begin


        test_P = momentum(psp, Incoming(), Electron())
        test_K = momentum(psp, Incoming(), Photon())

        test_P_prime = momentum(psp, Outgoing(), Electron())
        test_K_prime = momentum(psp, Outgoing(), Photon())

        @test isapprox(getMass2(test_P), mass(Electron())^2, atol = ATOL, rtol = RTOL)
        @test isapprox(getMass2(test_K), mass(Photon())^2, atol = ATOL, rtol = RTOL)
        @test isapprox(getMass2(test_P_prime), mass(Electron())^2, atol = ATOL, rtol = RTOL)
        @test isapprox(getMass2(test_K_prime), mass(Photon())^2, atol = ATOL, rtol = RTOL)
    end

    @testset "differential cross section" begin
        test_diffCS = differential_cross_section(psp)
        groundtruth = _groundtruth_diffCS(cth)

        @test isapprox(test_diffCS, groundtruth)
    end
end
end # ThomsonTests
