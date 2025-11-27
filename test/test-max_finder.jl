module MaxFinderTest

using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using XRTSProbing

#RNG = Xoshiro(137)
RNG = Random.Xoshiro(0xade95f7911a5ed27, 0x6afdf5bc148e56ef, 0x889de92d3feffd84, 0x3434b2628160aa75, 0x0769af4118a58ebb)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

const PROC = Thomson(PolX(), PolX())
const MODEL = PerturbativeQED()
const INPSL = TwoBodyTargetSystem()
const OUTPSL = PhotonSphericalLayout(INPSL)

const NCALLS = 50000
const MAXITER = 50

const QUANTILES = (0.01, 0.001, 0.0001)

const OMEGAS = (rand(RNG), 1.0e2 * rand(RNG), 1.0e3 * rand(RNG), 1.0e4 * rand(RNG))


@testset "om: $om" for om in OMEGAS
    DCSCACHED = DifferentialCrossSection(PROC, MODEL, OUTPSL; energy_2 = om)


    # build proposal
    VP = VegasProposal(DCSCACHED)

    train!(RNG, VP, MAXITER, NCALLS)

    @testset "p: $p" for p in QUANTILES
        test_max_finder = QuantileReductionMethod(p, Int(1.0e6))

        test_max_weight = XRTSProbing.findmax(RNG, DCSCACHED, test_max_finder, VP)

        @testset "n: $n" for n in (Int(1.0e6),)
            # groundtruth
            samples, jac = XRTSProbing._generate_coords(RNG, VP, n)
            weights = sort(@. DCSCACHED(samples) * jac)
            resid_weights = @. max(1, weights / test_max_weight)
            idx_last_unit_weight =
                length(resid_weights) - length(resid_weights[resid_weights .> 1.0])

            test_quantile =
                1.0 - sum(resid_weights[1:idx_last_unit_weight]) / sum(resid_weights)

            @test p >= test_quantile
        end
    end
end
end
