using Test
using SafeTestsets

begin
    @safetestset "Thomson process" begin
        include("hard_scattering/thomson/process.jl")
    end

    @safetestset "Thomson kinematics" begin
        include("hard_scattering/thomson/kinematics.jl")
    end

    @safetestset "Thomson cross section" begin
        include("hard_scattering/thomson/cross_section.jl")
    end
end
