using Test
using Random

using QEDcore
using QEDbase
using QEDprocesses
using QEDprobing

RNG = Xoshiro(137)
ATOL = 2 * sqrt(eps())
RTOL = sqrt(eps())

include("groundtruths.jl")
include("../checks.jl")

const Es = (1.0, 1.0 + rand(RNG), 2.0 + rand(RNG), 10.0 + rand(RNG))
const OMEGAS = (rand(RNG), 1.0e2 * rand(RNG), rand(RNG), 1.0e3 * rand(RNG), 1.0e4 * rand(RNG))
const CTHS = (-1.0, 0.0, 1.0, rand(RNG), -rand(RNG))
const PHIS = [0.0, 2.0 * pi, rand(RNG) * 2 * pi]

const MODEL = PerturbativeQED()
PROC = Thomson()
INPSLS = (TwoBodyTargetSystem(), PhotonElectronHeadsOnSystem())
#INPSL = TwoBodyTargetSystem()

KINPARAMETER = Dict(
    TwoBodyTargetSystem() => (
        (1.0,),
        OMEGAS,
        CTHS,
        PHIS,
    ),
    PhotonElectronHeadsOnSystem() => (
        Es,
        OMEGAS,
        CTHS,
        PHIS,
    )
)


@testset "kinmode: $kmode" for kmode in (Elastic(), InElastic())

    @testset "in_psl = $in_psl" for in_psl in INPSLS

        OUTPSL = PhotonSphericalLayout(in_psl, kmode)
        @testset "($E,$om,$cth,$phi)" for (E, om, cth, phi) in
            Iterators.product(KINPARAMETER[in_psl]...)

            # TODO: improve this!
            in_ps = in_psl == TwoBodyTargetSystem() ? (om,) : (E, om)
            out_ps = (cth, phi)
            psp = PhaseSpacePoint(PROC, MODEL, OUTPSL, in_ps, out_ps)

            @testset "particle momenta" begin

                test_P = momentum(psp, Incoming(), Electron())
                test_K = momentum(psp, Incoming(), Photon())

                test_P_prime = momentum(psp, Outgoing(), Electron())
                test_K_prime = momentum(psp, Outgoing(), Photon())

                @testset "onshell check" begin
                    _onshell_check(Electron(), test_P, ATOL, RTOL)
                    _onshell_check(Photon(), test_K, ATOL, RTOL)

                    if !is_elastic(OUTPSL)
                        _onshell_check(Electron(), test_P_prime, ATOL, RTOL)
                        _onshell_check(Photon(), test_K_prime, ATOL, RTOL)
                    end

                end

                @testset "components" begin
                    if cth^2 != one(cth)
                        @test isapprox(getE(test_K), om)
                        @test isapprox(getCosPhi(test_K_prime), cos(phi))
                        @test isapprox(getCosTheta(test_K_prime), cth)
                    end
                end


            end
        end
    end
end
