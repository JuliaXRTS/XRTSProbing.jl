using Test
using SafeTestsets

begin
    @safetestset "Lindhard" begin
        include("structure_factor/dynamic/lindhard.jl")
    end
end
