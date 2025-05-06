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

const MODEL = PerturbativeQED()
const INPSL = TwoBodyTargetSystem()

@testset "static properties" begin
    @testset "unpolarized" begin

        @test Thomson() == Thomson(AllSpin(), AllPol(), AllSpin(), AllPol())

        @testset "$INPOL" for INPOL in (PolX(), PolY())

            @test Thomson(INPOL) == Thomson(AllSpin(), INPOL, AllSpin(), AllPol())

            @testset "$OUTPOL" for OUTPOL in (PolX(), PolY())
                @test Thomson(INPOL, OUTPOL) == Thomson(AllSpin(), INPOL, AllSpin(), OUTPOL)
            end

        end
    end

    @testset "polarized" begin
        @testset "$INPOL, $OUTPOL" for (INPOL, OUTPOL) in (
            (AllPol(), AllPol()),
            Iterators.product((PolX(), PolY()), (PolX(), PolY()))...,
        )
            PROC = Thomson(INPOL, OUTPOL)

            @test incoming_particles(PROC) == (Electron(), Photon())
            @test outgoing_particles(PROC) == (Electron(), Photon())

            @test QEDbase.in_phase_space_dimension(PROC, MODEL) == 1
            @test QEDbase.out_phase_space_dimension(PROC, MODEL) == 2

        end
    end

    @testset "kinematic modes" begin
        @test is_elastic(Elastic()) == true
        @test is_elastic(InElastic()) == false

        test_psl_elastic = PhotonSphericalLayout(INPSL, Elastic())
        @test is_elastic(test_psl_elastic) == true

        test_psl_inelastic = PhotonSphericalLayout(INPSL, InElastic())
        @test is_elastic(test_psl_inelastic) == false

    end
end

@testset "kinmode: $kmode" for kmode in (Elastic(), InElastic())
    @testset "$INPOL, $OUTPOL" for (INPOL, OUTPOL) in (
        (AllPol(), AllPol()),
        Iterators.product((PolX(), PolY()), (PolX(), PolY()))...,
    )

        PROC = Thomson(INPOL, OUTPOL)
        OUTPSL = PhotonSphericalLayout(INPSL, kmode)

        @testset "($om,$cth,$phi)" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, CTHS, PHIS)
            in_ps = (om,)
            out_ps = (cth, phi)
            psp = PhaseSpacePoint(PROC, MODEL, OUTPSL, in_ps, out_ps)

            @testset "particle momenta" begin

                test_P = momentum(psp, Incoming(), Electron())
                test_K = momentum(psp, Incoming(), Photon())

                test_P_prime = momentum(psp, Outgoing(), Electron())
                test_K_prime = momentum(psp, Outgoing(), Photon())

                @test isapprox(
                    getMass2(test_P),
                    mass(Electron())^2,
                    atol = ATOL,
                    rtol = RTOL,
                )
                @test isapprox(getMass2(test_K), mass(Photon())^2, atol = ATOL, rtol = RTOL)

                if !is_elastic(OUTPSL)
                    @test isapprox(
                        getMass2(test_P_prime),
                        mass(Electron())^2,
                        atol = ATOL,
                        rtol = RTOL,
                    )
                    @test isapprox(
                        getMass2(test_K_prime),
                        mass(Photon())^2,
                        atol = ATOL,
                        rtol = RTOL,
                    )
                end

                @testset "components" begin
                    if cth^2 != one(cth)
                        @test isapprox(getE(test_K), om)
                        @test isapprox(getCosPhi(test_K_prime), cos(phi))
                        @test isapprox(getCosTheta(test_K_prime), cth)
                    end
                end
            end

            @testset "differential cross section" begin
                groundtruth = _groundtruth_TS_diffCS(INPOL, OUTPOL, cth, phi)

                @testset "direct" begin
                    test_diffCS = differential_cross_section(psp)
                    @test isapprox(test_diffCS, groundtruth)
                end

                @testset "from mat_el squared" begin
                    test_diffCS_from_MES = prod((
                        inv(4 * QEDbase._incident_flux(psp)),
                        QEDbase._matrix_element_square(psp),
                        QEDbase._averaging_norm(Float64, process(psp)),
                        QEDbase._phase_space_factor(psp),
                    ))
                    @test isapprox(test_diffCS_from_MES, groundtruth)
                end
            end
        end
    end
end
end # ThomsonTests
