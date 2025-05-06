using Test
using Random
using Unitful

using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

NE = rand(RNG) * 1.0e21u"cm^(-3)"
#NE = 1e21u"cm^(-3)"
NE_internal = QEDprobing._internalize_density(NE)
KF = QEDprobing._fermi_wave_vector(NE_internal)
EF = QEDprobing._fermi_energy_from_kF(KF)

OMS = EF .* (0.0, rand(RNG), 1 + rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
QS = KF .* (1.0e-2 * rand(RNG), rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
TEMPS_eV = (rand(RNG), rand(RNG) * 10, rand(RNG) * 100) .* 1u"eV"
TEMPS = QEDprobing._internalize_temperature.(TEMPS_eV)

TEST_SYSTEM_zeroT = IdealElectronSystem{ZeroTemperature}(NE)
TEST_SYSTEMS_finT = IdealElectronSystem.(NE, TEMPS)
TEST_SYSTEMS = (TEST_SYSTEM_zeroT, TEST_SYSTEMS_finT...)


@testset "temp: $T" for (i, T) in enumerate((zero(eltype(TEMPS)), TEMPS...))
    test_system = TEST_SYSTEMS[i]

    @testset "no screening" begin

        scr = NoScreening()
        @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)

            @test isapprox(dielectric_function(test_system, scr, (om, q)), one(om))

            @test isapprox(local_effective_potential(scr, (om, q)), zero(om))

            @test isapprox(pseudo_potential(scr, (om, q)), zero(om))

            @test isapprox(local_field_correction(scr, (om, q)), zero(om))
        end
    end

    @testset "screening" begin

        # checking only default value: coulomb+no_lfc
        scr = Screening()

        @testset "default constructor" begin
            @test scr == Screening(CoulombPseudoPotential(), NoLocalFieldCorrection())
        end

        @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)

            @test isapprox(
                dielectric_function(test_system, scr, (om, q)),
                one(om) -
                    pseudo_potential(scr, (om, q)) * dynamic_response(test_system, (om, q)),
            )

            @test isapprox(
                local_effective_potential(scr, (om, q)),
                pseudo_potential(scr, (om, q)),
            )

            @test isapprox(pseudo_potential(scr, (om, q)), ELEMENTARY_CHARGE_SQUARED / q^2)

            @test isapprox(local_field_correction(scr, (om, q)), zero(om))
        end
    end
end
