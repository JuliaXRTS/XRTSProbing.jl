using Test
using SafeTestsets

begin
    @safetestset "abstract setup" begin
        include("setups/abstract.jl")
    end

    @safetestset "differential cross section" begin
        include("setups/diff_cross_section.jl")
    end

    @safetestset "probing setup" begin
        include("setups/probing.jl")
    end
end
