using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using XRTSProbing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

include("groundtruths.jl")

const MODEL = PerturbativeQED()
#INPSLS = (TwoBodyTargetSystem(),PhotonElectronHeadsOnSystem())
INPSL = TwoBodyTargetSystem()


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
