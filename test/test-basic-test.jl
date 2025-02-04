module QEDprobingTest

using Test
using QEDprobing


@testset "QEDprobing.jl" begin
    @test QEDprobing.hello_world() == "Hello, World!"
end

end
