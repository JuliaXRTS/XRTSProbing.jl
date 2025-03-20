using QEDprobing
using QEDcore
using QEDevents

using Distributions
using Random
RNG = Xoshiro(137)

include("utils.jl")

const AXES = (XAxis(), YAxis(), ZAxis(), AllAxis())
const OMEGA_MEANS = (
    1e-4 * rand(RNG),
    1e-3 * rand(RNG),
    1e-2 * rand(RNG),
    1e-1 * rand(RNG),
    rand(RNG),
    1e2 * rand(RNG),
    rand(RNG),
    1e3 * rand(RNG),
    1e4 * rand(RNG),
)
const OMEGA_WIDTHS = (
    0.0,
    1e-2 * rand(RNG),
    1e-3 * rand(RNG),
    1e-4 * rand(RNG),
    1e-1 * rand(RNG),
    rand(RNG),
    1e2 * rand(RNG),
    rand(RNG),
    1e3 * rand(RNG),
    1e4 * rand(RNG),
)

@testset "om: $om, width: $width" for (om, width) in
                                      Iterators.product(OMEGA_MEANS, OMEGA_WIDTHS)

    @testset "default" begin
        DIST = GaussianPhotonDist(om, width)
        @test DIST == GaussianPhotonDist(om, width, AllAxis())
    end

    @testset "$axis" for axis in AXES
        DIST = GaussianPhotonDist(om, width, axis)
        groundtruth_energy_dist = truncated(Normal(om, width), 0.0, Inf)

        @testset "properties" begin
            @test QEDevents._particle(DIST) == Photon()
            @test QEDevents._particle_direction(DIST) == Incoming()
            @test energy_mean(DIST) == om
            @test energy_width(DIST) == width
            @test k_vec_axis(DIST) == axis
        end

        @testset "weights" begin
            rnd_omega = rand(RNG, groundtruth_energy_dist)
            groundtruth_weight = pdf(groundtruth_energy_dist, rnd_omega)
            @test isapprox(weight(DIST, rnd_omega), groundtruth_weight)
            @test isapprox(weight(DIST, -rnd_omega), 0.0)
        end

        @testset "random momenta" begin

            # TODO: test `randmom` here if available
            rnd_mom = QEDevents._randmom(RNG, DIST)
            if axis isa QEDprobing.AbstractDefiniteAxis
                @test isapprox(getMass2(rnd_mom), 0.0)

                for (i, el) in enumerate(rnd_mom)
                    if i == _axis_index(axis) || i == 1
                        @test el != 0.0
                    else
                        @test isapprox(el, 0.0)
                    end
                end

            else # AllAxis case
                # TODO: figure out, why this case has some precision problems.
                @test _is_onshell_photon(rnd_mom)
            end
        end

    end

end
