
module SetupTest

using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())


function _get_all_except(t::NTuple{N,T}, exc::T) where {N,T}
    t_a = [t...]
    return Tuple(t_a[t_a.!=exc])
end

const OMEGAS = (rand(RNG), 1e2 * rand(RNG), rand(RNG), 1e3 * rand(RNG), 1e4 * rand(RNG))
const CTHS = (-1.0, 0.0, 1.0, rand(RNG), -rand(RNG))
const PHIS = (0.0, 2 * pi, rand(RNG) * 2 * pi)
const COORD_SYMS = (:energy_2, :cos_theta, :phi)
const COORDS = (OMEGAS, CTHS, PHIS)

const MODEL = PerturbativeQED()
const INPSL = TwoBodyTargetSystem()
const OUTPSL = PhotonSphericalLayout(INPSL)

@testset "$INPOL, $OUTPOL" for (INPOL, OUTPOL) in (
    (AllPol(), AllPol()),
    Iterators.product((PolX(), PolY()), (PolX(), PolY()))...,
)
    PROC = Thomson(INPOL, OUTPOL)

    DIFFCS_SETUP = DifferentialCrossSection(PROC, MODEL, OUTPSL)
    @testset "full input" begin


        @testset "($om,$cth,$phi)" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, CTHS, PHIS)
            in_ps = (om,)
            out_ps = (cth, phi)
            all_ps = (om, cth, phi)
            psp = PhaseSpacePoint(PROC, MODEL, OUTPSL, in_ps, out_ps)

            groundtruth_diff_cs = differential_cross_section(psp)

            diff_cs = DIFFCS_SETUP(in_ps, out_ps)
            diff_cs_psp = DIFFCS_SETUP(psp)
            diff_cs_all_ps = DIFFCS_SETUP(all_ps)

            @test isapprox(diff_cs, groundtruth_diff_cs, atol = ATOL, rtol = RTOL)
            @test isapprox(diff_cs_psp, groundtruth_diff_cs, atol = ATOL, rtol = RTOL)
            @test isapprox(diff_cs_all_ps, groundtruth_diff_cs, atol = ATOL, rtol = RTOL)
        end
    end

    @testset "kw cached input" begin
        @testset "single kw" begin
            @testset "$sym" for (i, sym) in enumerate(COORD_SYMS)

                @testset "($om,$cth,$phi)" for (om, cth, phi) in
                                               Iterators.product(OMEGAS, CTHS, PHIS)
                    test_coords = [om, cth, phi]
                    test_cached_coord = test_coords[i]
                    test_rem_coords = Tuple(test_coords[1:end.!=i])

                    test_diff_cs_setup = DifferentialCrossSection(
                        PROC,
                        MODEL,
                        OUTPSL;
                        sym => test_cached_coord,
                    )
                    test_diff_cs = test_diff_cs_setup(test_rem_coords)
                    groundtruth = DIFFCS_SETUP(Tuple(test_coords))

                    @test isapprox(test_diff_cs, groundtruth)
                end
            end
        end
        @testset "multiple kw" begin
            for (i1, sym1) in enumerate(COORD_SYMS)
                for (i2, sym2) in enumerate(COORD_SYMS)
                    if sym1 != sym2
                        @testset "$sym1 $sym2" begin
                            @testset "($om,$cth,$phi)" for (om, cth, phi) in
                                                           Iterators.product(
                                OMEGAS,
                                CTHS,
                                PHIS,
                            )
                                test_coords = [om, cth, phi]
                                test_cached_coord =
                                    (; sym1 => test_coords[i1], sym2 => test_coords[i2])
                                test_rem_coords =
                                    Tuple(test_coords[(1:end.!=i1).&&(1:end.!=i2)])

                                test_diff_cs_setup = DifferentialCrossSection(
                                    PROC,
                                    MODEL,
                                    OUTPSL;
                                    test_cached_coord...,
                                )
                                test_diff_cs = test_diff_cs_setup(test_rem_coords)
                                groundtruth = DIFFCS_SETUP(Tuple(test_coords))

                                @test isapprox(test_diff_cs, groundtruth)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "cached input" begin
        @testset "($om,$cth,$phi)" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, CTHS, PHIS)
            in_ps = (om,)
            out_ps = (cth, phi)

            DIFFCS_SETUP = DifferentialCrossSectionCached(PROC, MODEL, OUTPSL, in_ps)

            psp = PhaseSpacePoint(PROC, MODEL, OUTPSL, in_ps, out_ps)

            groundtruth_diff_cs = differential_cross_section(psp)

            diff_cs = DIFFCS_SETUP(out_ps)

            @test isapprox(diff_cs, groundtruth_diff_cs, atol = ATOL, rtol = RTOL)
        end
    end
end
end
