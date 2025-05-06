using QEDprobing
using QEDcore
using QEDevents

using Distributions
using Random
RNG = Xoshiro(137)

include("utils.jl")

const AXES = (XAxis(), YAxis(), ZAxis(), UniformAxis())
const OMEGA_MEANS = (
    1.0e-4 * rand(RNG),
    1.0e-3 * rand(RNG),
    1.0e-2 * rand(RNG),
    1.0e-1 * rand(RNG),
    rand(RNG),
    1.0e2 * rand(RNG),
    rand(RNG),
    1.0e3 * rand(RNG),
    1.0e4 * rand(RNG),
)
const OMEGA_WIDTHS = (
    0.0,
    1.0e-2 * rand(RNG),
    1.0e-3 * rand(RNG),
    1.0e-4 * rand(RNG),
    1.0e-1 * rand(RNG),
    rand(RNG),
    1.0e2 * rand(RNG),
    rand(RNG),
    1.0e3 * rand(RNG),
    1.0e4 * rand(RNG),
)

@testset "$axis" for axis in AXES
    @testset "om: $om, width: $width" for (om, width) in
        Iterators.product(OMEGA_MEANS, OMEGA_WIDTHS)
        photon_dist = GaussianPhotonDist(om, width, axis)

        photon_field = DistributionBasedField(photon_dist)

        test_omegas = rand(RNG, photon_dist.dist, 10)

        @test isapprox(
            energy_spectrum.(photon_field, test_omegas),
            weight.(photon_dist, test_omegas),
        )

    end

end
