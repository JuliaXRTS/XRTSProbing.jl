module RejectionTest

using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

const PROC = Thomson(PolX(), PolX())
const MODEL = PerturbativeQED()
const INPSL = TwoBodyTargetSystem()
const OUTPSL = PhotonSphericalLayout(INPSL)

const NCALLS = 50000
const MAXITER = 50
const MAXFINDER = QuantileReductionMethod(0.001, Int(1.0e6))

const NEVENTS = (1, rand(RNG, 2:10), rand(RNG, 11:20))
const BATCHSIZE = (1, rand(RNG, 2:10), rand(RNG, 11:20))

const OMEGAS = (rand(RNG), 1.0e2 * rand(RNG), rand(RNG), 1.0e3 * rand(RNG), 1.0e4 * rand(RNG))


@testset "om: $om" for om in OMEGAS
    test_in_coords = (om,)
    DCSCACHED = DifferentialCrossSectionCached(PROC, MODEL, OUTPSL, test_in_coords)


    # build proposal
    VP = VegasProposal(DCSCACHED)

    train!(RNG, VP, MAXITER, NCALLS)

    # max finding
    test_max_weight = QEDprobing.findmax(RNG, DCSCACHED, MAXFINDER, VP)

    # build event generator
    EG = EventGenerator(DCSCACHED, VP, test_max_weight)

    @testset "properties" begin
        @test process(EG) == PROC
        @test model(EG) == MODEL
        @test phase_space_layout(EG) == OUTPSL
    end

    @testset "generation" begin

        @testset "single" begin
            event = generate_event(RNG, EG)

            @test process(event.psp) == PROC
            @test model(event.psp) == MODEL
            @test phase_space_layout(event.psp) == OUTPSL

            @test event.weight >= one(event.weight)
        end

        @testset "multiple" begin
            @testset "N: $n, bs: $bs" for (n, bs) in Iterators.product(NEVENTS, BATCHSIZE)
                event_list = generate_events(RNG, EG, n, bs)

                @test length(event_list) == n
                for event in event_list
                    @test process(event.psp) == PROC
                    @test model(event.psp) == MODEL
                    @test phase_space_layout(event.psp) == OUTPSL

                    @test event.weight >= one(event.weight)
                end
            end
        end


    end


end
end
