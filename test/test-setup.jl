
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


const OMEGAS = (rand(RNG), 1e2 * rand(RNG), rand(RNG), 1e3 * rand(RNG), 1e4 * rand(RNG))
const CTHS = (-1.0, 0.0, 1.0, rand(RNG), -rand(RNG))
const PHIS = [0, 2 * pi, rand(RNG) * 2 * pi]

const MODEL = PerturbativeQED()
const INPSL = TwoBodyTargetSystem()
const OUTPSL = PhotonSphericalLayout(INPSL)

@testset "$INPOL, $OUTPOL" for (INPOL, OUTPOL) in (
    (AllPol(), AllPol()),
    Iterators.product((PolX(), PolY()), (PolX(), PolY()))...,
)
    PROC = Thomson(INPOL, OUTPOL)

    @testset "full input" begin

        DIFFCS_SETUP = DifferentialCrossSection(PROC, MODEL, OUTPSL)

        @testset "($om,$cth,$phi)" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, CTHS, PHIS)
            in_ps = (om,)
            out_ps = (cth, phi)
            psp = PhaseSpacePoint(PROC, MODEL, OUTPSL, in_ps, out_ps)

            groundtruth_diff_cs = differential_cross_section(psp)

            diff_cs = DIFFCS_SETUP(in_ps, out_ps)

            @test isapprox(diff_cs, groundtruth_diff_cs, atol = ATOL, rtol = RTOL)
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
