
module KinCutsTests

using Test
using QEDprobing
using Random
using Distributions

RNG = Xoshiro(137)

DOFs = [1, rand(RNG, 2:10)]

@testset "dof: $dof" for dof in DOFs
    low = Tuple(-rand(RNG, dof))
    up = Tuple(rand(RNG, dof))

    cuts = KinematicCuts(low, up)

    @testset "properties" begin
        @test degree_of_freedom(cuts) == dof
        @test minimum(cuts) == low
        @test maximum(cuts) == up
        @test extrema(cuts) == (low, up)
    end

    @testset "element check" begin
        all_valid_dist = product_distribution(Uniform.(low, up)...)
        all_valid_element = Tuple(rand(RNG, all_valid_dist))

        @test all_within_cuts(cuts, all_valid_element)
        @test all(is_within_cuts(cuts, all_valid_element))

        all_invalid_dist = product_distribution(Uniform.(2 .* low, low)...)
        all_invalid_element = Tuple(rand(RNG, all_invalid_dist))

        @test all_within_cuts(cuts, all_invalid_element) == false
        @test all(is_within_cuts(cuts, all_invalid_element)) == false

        if dof > 1
            # switching one element in validity makes no sense of dof==1

            one_valid_element = (0.0, all_invalid_element[2:end]...)
            groundtruth_one_valid = tuple(true, fill(false, dof - 1)...)
            @test all_within_cuts(cuts, one_valid_element) == false
            @test is_within_cuts(cuts, one_valid_element) == groundtruth_one_valid

            one_invalid_element = (2 * low[1], all_valid_element[2:end]...)
            groundtruth_one_invalid = .!groundtruth_one_valid
            @test all_within_cuts(cuts, one_valid_element) == false
            @test is_within_cuts(cuts, one_invalid_element) == groundtruth_one_invalid
        end
    end

    @testset "construction" begin
        # exchange all bounds
        @test_throws ArgumentError KinematicCuts(up, low)

        # exchange one last bounds
        low_one_invalid = (low[1:end-1]..., up[end])
        up_one_invalid = (up[1:end-1]..., low[end])
        @test_throws ArgumentError KinematicCuts(low_one_invalid, up_one_invalid)

    end

end

end
