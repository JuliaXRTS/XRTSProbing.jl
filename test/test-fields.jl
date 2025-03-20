
using Test
using SafeTestsets

begin
    @safetestset "photon dists" begin
        include("fields/photon_dists.jl")
    end

    @safetestset "dist based field" begin
        include("fields/dist_based.jl")
    end
end
