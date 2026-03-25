using Test
using SafeTestsets

begin
    @safetestset "Maxwell Boltzmann" begin
        include("electron_dists/maxwell_dist.jl")
    end

    @safetestset "Uniform" begin
        include("electron_dists/uniform_dist.jl")
    end
end
