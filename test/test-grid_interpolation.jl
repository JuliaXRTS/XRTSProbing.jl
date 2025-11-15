using QEDprobing
using Random
RNG = Xoshiro(137)
gaussian(x, y) = exp(-x^2 / 2 - y^2)
a_arr = collect(range(-rand(RNG), rand(RNG); length = 4)) # cols
b_arr = collect(range(-rand(RNG), rand(RNG); length = 5)) # rows
vals = gaussian.(a_arr', b_arr)
interpol_methods = [
    InterpolExtrapol(),
    InterpolEndValue(),
    NearestInput(),
    BelowInput(),
    AboveInput(),
]


@testset "$method" for method in interpol_methods

    @testset "sanity check" begin
        @test GridInterpolant(a_arr, b_arr, vals, method) == interpolate(a_arr, b_arr, vals, method)
    end

    @testset "reproduction" begin
        ITP = interpolate(a_arr, b_arr, vals, method)
        for a in a_arr
            for b in b_arr
                test_val = ITP(a, b)
                groundtruth = gaussian(a, b)
                @test isapprox(groundtruth, test_val)
            end
        end
    end

end
