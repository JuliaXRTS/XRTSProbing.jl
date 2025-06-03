module MaxFinderTest

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

const QUANTILES = (0.01, 0.001, 0.0001)

const OMEGAS = (rand(RNG), 1.0e2 * rand(RNG), rand(RNG), 1.0e3 * rand(RNG), 1.0e4 * rand(RNG))


@testset "om: $om" for om in OMEGAS
    DCSCACHED = DifferentialCrossSection(PROC, MODEL, OUTPSL; energy_2 = om)


    # build proposal
    VP = VegasProposal(DCSCACHED)

    train!(RNG, VP, MAXITER, NCALLS)

    @testset "p: %p" for p in QUANTILES
        test_max_finder = QuantileReductionMethod(p, Int(1.0e6))

        test_max_weight = QEDprobing.findmax(RNG, DCSCACHED, test_max_finder, VP)

        @testset "n: %n" for n in (Int(2.0e5), Int(1.0e6))
            # groundtruth
            samples, jac = QEDprobing._generate_coords(RNG, VP, n)
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
