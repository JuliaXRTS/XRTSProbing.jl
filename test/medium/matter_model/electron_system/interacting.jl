using Test
using Random
using QuadGK
using Unitful

using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = 1.0e-6

NE = rand(RNG) * 1.0e21u"cm^(-3)"
#NE = 1e21u"cm^(-3)"
NE_internal = QEDprobing._internalize_density(NE)
KF = QEDprobing._fermi_wave_vector(NE_internal)
EF = QEDprobing._fermi_energy_from_kF(KF)

OMS = EF .* (0.0, rand(RNG), 1 + rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
QS = KF .* (1.0e-2 * rand(RNG), rand(RNG), 2 + rand(RNG), 3 + rand(RNG))
TEMPS_eV = (rand(RNG), rand(RNG) * 10, rand(RNG) * 100) .* 1u"eV"
TEMPS = QEDprobing._internalize_temperature.(TEMPS_eV)

TEST_PROPER_SYSTEM_zeroT = IdealElectronSystem{ZeroTemperature}(NE)
TEST_PROPER_SYSTEMS_finT = IdealElectronSystem.(NE, TEMPS)
TEST_PROPER_SYSTEMS = (TEST_PROPER_SYSTEM_zeroT, TEST_PROPER_SYSTEMS_finT...)

SCREENINGS = (Screening(), NoScreening())

@testset "$scr" for scr in SCREENINGS

    @testset "zero temperature" begin
        test_system_zeroT = InteractingElectronSystem(TEST_PROPER_SYSTEM_zeroT, scr)
        @testset "properties" begin
            @test iszero(temperature(test_system_zeroT))
            @test isinf(beta(test_system_zeroT))
            @test isinf(betabar(test_system_zeroT))
        end

        @testset "matching" begin
            test_system_small_finT =
                InteractingElectronSystem(IdealElectronSystem(NE, 1.0e-2 * eps()), scr)

            @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)
                @test isapprox(
                    imag_dynamic_response(test_system_zeroT, (om, q)),
                    imag_dynamic_response(test_system_small_finT, (om, q)),
                    rtol = 1.0e-2,
                    atol = 1.0e-6,
                )
            end
        end
    end

    @testset "temp: $T" for (i, T) in enumerate((zero(eltype(TEMPS)), TEMPS...))
        test_system = InteractingElectronSystem(TEST_PROPER_SYSTEMS[i], scr)

        @testset "properties" begin
            @test temperature(test_system) == T
            @test isapprox(beta(test_system), inv(T))
            @test fermi_wave_vector(test_system) == KF
            @test fermi_energy(test_system) == EF
            @test electron_density(test_system) == NE_internal
        end

        @testset "stability" begin
            @testset "q: $q" for q in QS
                @test real_dynamic_response(test_system, (0.0, q)) <= zero(q)
                @test imag_dynamic_response(test_system, (0.0, q)) <= zero(q)
            end
        end

        @testset "sanity check" begin
            @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)

                groundtruth_imag_rf = imag(dynamic_response(test_system, (om, q)))
                @test isapprox(
                    groundtruth_imag_rf,
                    imag_dynamic_response(test_system, (om, q)),
                    rtol = RTOL,
                )

                groundtruth_real_rf = real(dynamic_response(test_system, (om, q)))

                @test isapprox(
                    groundtruth_real_rf,
                    real_dynamic_response(test_system, (om, q)),
                    rtol = RTOL,
                )

                groundtruth_rf =
                    real_dynamic_response(test_system, (om, q)) +
                    1im * imag_dynamic_response(test_system, (om, q))
                @test isapprox(
                    groundtruth_rf,
                    dynamic_response(test_system, (om, q)),
                    rtol = RTOL,
                )

                proper_rf = dynamic_response(TEST_PROPER_SYSTEMS[i], (om, q))
                lep = local_effective_potential(test_system, (om, q))
                groundtruth_rf_from_proper = proper_rf / (1 - lep * proper_rf)
                @test isapprox(
                    groundtruth_rf_from_proper,
                    dynamic_response(test_system, (om, q)),
                    rtol = RTOL,
                )

                diel_func = dielectric_function(test_system, (om, q))
                groundtruth_rf_from_dielec_func = proper_rf / diel_func
                @test isapprox(
                    groundtruth_rf_from_dielec_func,
                    dynamic_response(test_system, (om, q)),
                    rtol = RTOL,
                )
            end
        end

        @testset "symmetry" begin
            @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)
                input_pos_q = (om, q)
                input_pos_om = (om, q)
                input_neg_om = (-om, q)
                input_neg_both = (-om, -q)
                @test isapprox(
                    dynamic_response(test_system, input_pos_q),
                    dynamic_response(test_system, input_neg_both),
                    rtol = RTOL,
                )
                @test isapprox(
                    real_dynamic_response(test_system, input_pos_om),
                    real_dynamic_response(test_system, input_neg_om),
                    rtol = RTOL,
                )
                @test isapprox(
                    imag_dynamic_response(test_system, input_pos_om),
                    -imag_dynamic_response(test_system, input_neg_om),
                    rtol = RTOL,
                )
            end
        end
        #=
            @testset "f-sum rule" begin
                if iszero(T)
                    @testset "q: $q" for q in QS
                        qb = q/KF
                        lower_om_bound = EF*max(zero(qb),qb^2 - 2 * qb)
                        upper_om_bound = EF*(qb^2 + 2 * qb)
                        tmp, err = quadgk(
                            x -> x * imag_dynamic_response(test_system,(x, q)),
                            lower_om_bound,
                            upper_om_bound
                        )
                        _first_moment = -2 * tmp / pi
                        @test isapprox(_first_moment, 2*q^2*fermi_energy(test_system)/3,rtol = 1e-2)
                    end
                end
            end
                =#

        @testset "dynamic structure factor" begin
            @testset "om: $om, q: $q" for (om, q) in Iterators.product(OMS, QS)
                @testset "sanity check" begin
                    if !iszero(om)
                        fac = pi * electron_density(test_system)

                        im_rf = imag_dynamic_response(test_system, (om, q))
                        groundtruth_dsf =
                            iszero(T) ? -im_rf / fac :
                            -inv(fac * (one(om) - exp(-om / temperature(test_system)))) *
                            im_rf

                        dsf = dynamic_structure_factor(test_system, (om, q))
                        @test isapprox(dsf, groundtruth_dsf)
                        @test dsf >= zero(dsf)
                    end
                end

                @testset "detailed balance" begin
                    if !iszero(om) && !iszero(T)
                        input_pos_om = (om, q)
                        input_neg_om = (-om, q)

                        @test isapprox(
                            dynamic_structure_factor(test_system, input_neg_om),
                            exp(-om * beta(test_system)) *
                                dynamic_structure_factor(test_system, input_pos_om),
                            rtol = 1.0e-4,
                        )
                    end
                end
            end
        end
    end
end
