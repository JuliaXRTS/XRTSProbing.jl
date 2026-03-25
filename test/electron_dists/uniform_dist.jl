using Random
using Unitful

using XRTSProbing
using ElectronicStructureModels

RNG = Xoshiro(137)

DTYPES = (Float32, Float64)

@testset "dtype: $dtype" for dtype in DTYPES
    test_dist = UniformElectronEnergyDistribution()

    en = rand(RNG, dtype)

    @test isapprox(XRTSProbing._energy_weight(test_dist, en), one(dtype))
end
