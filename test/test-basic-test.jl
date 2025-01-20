module QEDprobingTest

using Test
using QEDprobing


const A = 1

@testset "QEDprobing.jl" begin
    @test QEDprobing.hello_world() == "Hello, World!"
end

end
