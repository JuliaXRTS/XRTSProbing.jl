<<<<<<< HEAD
using Random
using Unitful

using QEDprobing

RNG = Xoshiro(137)

TEMPS_eV = 1u"eV" .* (
    rand(RNG),
    rand(RNG) * 10,
    rand(RNG) * 100,
    rand(RNG) * 1000,
)


@testset "temp = $T" for T in TEMPS_eV
    test_dist = MaxellElectronEnergyDistribution(T)

    T_internal = QEDprobing._internalize_temperature(T)

    @test isapprox(temperature(test_dist), T_internal)
    @test isapprox(energy_mean(test_dist), sqrt(8 * T_internal / pi + 1))
    @test MaxellElectronEnergyDistribution(T_internal) == test_dist

    # FIXME: enable if QEDevents.jl is fixed
    @test_broken isapprox(energy_width(test_dist), sqrt(T_internal * ((3 * pi - 8) / pi)^2 + 1))
end
=======
# TBW
>>>>>>> f4c8e3f (added electron dist; added event IO)
