using QEDprobing
using Test
using SafeTestsets

@testset "QEDprobing.jl" begin
    begin
        # patches
        @time @safetestset "patch QEDprocesses.jl" begin
            include("patch_QEDjl/patch_QEDprocesses.jl")
        end
    end
end
