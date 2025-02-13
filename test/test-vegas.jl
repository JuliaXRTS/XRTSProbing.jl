# implement some tests for vegas
#
# TODO: implement tests for
# - bin average tools (build,smooth, compress)
# - grid refinement
# - training procedure
# - sampling from the grid
#
using Test
using SafeTestsets

begin
    #=
    @safetestset "Grid" begin
        include("vegas/grid.jl")
    end

    @safetestset "Map" begin
        include("vegas/map.jl")
    end

    @safetestset "Jac" begin
        include("vegas/jac.jl")
    end
    =#
    @safetestset "bin avg" begin
        include("vegas/bin_avg.jl")
    end
end
