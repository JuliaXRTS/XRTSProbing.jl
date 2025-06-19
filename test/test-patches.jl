module Patches
using QEDprobing
using Test
using SafeTestsets

@safetestset "QEDevents.jl" begin
    include("patch_tests/QEDevents.jl")
end

@safetestset "QEDcore.jl" begin
    include("patch_tests/QEDcore.jl")
end

end
