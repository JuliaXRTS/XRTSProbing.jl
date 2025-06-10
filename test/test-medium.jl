using Test
using SafeTestsets

begin
    @safetestset "Ideal system" begin
        include("medium/electron_system/ideal.jl")
    end
    @safetestset "Screening" begin
        include("medium/electron_system/screening.jl")
    end
    @safetestset "Interacting system" begin
        include("medium/electron_system/interacting.jl")
    end
end
