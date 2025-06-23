using Test
using SafeTestsets

begin
    @safetestset "Thomson process" begin
        include("electron_dists/maxwell_dist.jl")
    end

end
