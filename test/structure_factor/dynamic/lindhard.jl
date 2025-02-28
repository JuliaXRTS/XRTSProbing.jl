using Test
using Random
using QuadGK

using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

OMBARS = (0.0, rand(RNG), 1 + rand(RNG), 10 + rand(RNG))
QBARS = (rand(RNG), 2 + rand(RNG), 10 + rand(RNG))
BETABARS = (rand(RNG), rand(RNG) * 10, rand(RNG) * 100)

TEST_RF_ZEROT = LindhardDRF()
TEST_RFS_FINT = LindhardDRF.(BETABARS)
TEST_RFS_ALL = (TEST_RFS_FINT..., TEST_RF_ZEROT)
@testset "response function" begin

    @testset "zero temperature" begin
        @test iszero(temperature(TEST_RF_ZEROT))
        @test isinf(betabar(TEST_RF_ZEROT))
    end

    @testset "finite temperature" begin
        @testset "betabar: $(betabar(test_rf))" for (idx, test_rf) in
                                                    enumerate(TEST_RFS_FINT)

            @test isapprox(betabar(test_rf), BETABARS[idx])
            @test isapprox(temperature(test_rf), inv(BETABARS[idx]))
            # TODO: add tests for other fields

        end
    end

    @testset "properties" begin
        @testset "betabar: $(betabar(test_rf))" for (idx, test_rf) in
                                                    enumerate(TEST_RFS_ALL)

            @testset "stability" begin
                @testset "qb: $qb" for qb in QBARS
                    @test real(test_rf, (0.0, qb)) <= zero(qb)
                    @test imag(test_rf, (0.0, qb)) <= zero(qb)
                end
            end

            @testset "symmetry" begin
                @testset "omb: $omb, qb: $qb" for (omb, qb) in
                                                  Iterators.product(OMBARS, QBARS)
                    input_pos_q = (omb, qb)
                    input_pos_om = (omb, qb)
                    input_neg_om = (-omb, qb)
                    input_neg_both = (-omb, -qb)
                    @test isapprox(
                        test_rf(input_pos_q),
                        test_rf(input_neg_both),
                        rtol = 1e-6,
                    )
                    @test isapprox(
                        real(test_rf, input_pos_om),
                        real(test_rf, input_neg_om),
                        rtol = 1e-6,
                    )
                    @test isapprox(
                        imag(test_rf, input_pos_om),
                        -imag(test_rf, input_neg_om),
                        rtol = 1e-6,
                    )
                end
            end

            @testset "f-sum rule" begin
                @testset "qb: $qb" for qb in QBARS

                    # bound only valid for positiv q, practical for 0.001<betabar<Inf
                    upper_omb_bound = qb^2 + 2 * qb + 80 / sqrt(betabar(test_rf))
                    tmp, err = quadgk(x -> x * imag(test_rf, (x, qb)), 0.0, upper_omb_bound)
                    _first_moment = -2 * tmp / pi
                    @test isapprox(_first_moment, 8 / 3 * qb^2, rtol = 1e-2)
                end
            end
        end
    end

    @testset "sanity check" begin

        # betabar=Inf corresponds to T=0
        test_rf_very_low_T = LindhardDRF(Inf)

        @testset "properties" begin
            @test isapprox(eps(), temperature(test_rf_very_low_T), atol = 0.0, rtol = 1e-6)
        end
        @testset "omb: $omb, qb: $qb" for (omb, qb) in Iterators.product(OMBARS, QBARS)

            # we allow worse approx, because the real part is an integral
            @test isapprox(
                test_rf_very_low_T((omb, qb)),
                TEST_RF_ZEROT((omb, qb)),
                rtol = 1e-4,
            )
        end
    end
end
