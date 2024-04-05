using QEDprobing
using Test

@testset "QEDprobing.jl" begin
    # Write your tests here.
    
    @testset "dummy" begin
        @test QEDprobing.dummy_func() == 1 
    end
end
