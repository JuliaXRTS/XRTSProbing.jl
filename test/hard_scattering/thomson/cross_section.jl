using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

include("groundtruths.jl")

const Es = (1.0,)
const OMEGAS = (rand(RNG), 1.0e2 * rand(RNG), rand(RNG), 1.0e3 * rand(RNG), 1.0e4 * rand(RNG))
const CTHS = (-1.0, 0.0, 1.0, rand(RNG), -rand(RNG))
const PHIS = [0.0, 2.0 * pi, rand(RNG) * 2 * pi]

const MODEL = PerturbativeQED()
#INPSLS = (TwoBodyTargetSystem(),PhotonElectronHeadsOnSystem())
INPSL = TwoBodyTargetSystem()

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

            @testset "differential cross section" begin
                groundtruth = _groundtruth_TS_diffCS(INPOL, OUTPOL, cth, phi)

                @testset "direct" begin
                    test_diffCS = differential_cross_section(psp)
                    @test isapprox(test_diffCS, groundtruth)
                end

                @testset "from mat_el squared" begin
                    test_diffCS_from_MES = prod(
                        (
                            inv(4 * QEDbase._incident_flux(psp)),
                            QEDbase._matrix_element_square(psp),
                            QEDbase._averaging_norm(Float64, process(psp)),
                            QEDbase._phase_space_factor(psp),
                        )
                    )
                    @test isapprox(test_diffCS_from_MES, groundtruth)
                end
            end
        end
    end
end
