using XRTSProbing
using QEDcore
using QEDprocesses
using Random

ATOL = 1.0e-12
RTOL = 1.0e-6

RNG = MersenneTwister(137137)
PARTICLES = (Electron(), Positron(), Photon())

N_INCOMING = 2
N_OUTGOING = rand(RNG, 2:8)

IN_PARTICLES = (Electron(), rand(RNG, PARTICLES))
OUT_PARTICLES = Tuple(rand(RNG, PARTICLES, N_INCOMING))

PROC = ScatteringProcess(IN_PARTICLES, OUT_PARTICLES)
MODEL = PertQED()
INPSL = TwoBodyRestSystem()
INCOORDS = (10.0,)

FLAT_PS_SAMPLER = FlatPhaseSpaceSampler(PROC, MODEL, INPSL, INCOORDS)
IN_MOMS, OUT_MOMS = randmom(RNG, FLAT_PS_SAMPLER)

@testset "energy momentum conservation" begin
    @test isapprox(sum(IN_MOMS), sum(OUT_MOMS))
end

@testset "on-shell check" begin
    for (mom, particle) in zip(IN_MOMS, IN_PARTICLES)
        @test isapprox(mom * mom, mass(particle), atol = ATOL, rtol = RTOL)
    end
    for (mom, particle) in zip(OUT_MOMS, OUT_PARTICLES)
        @test isapprox(mom * mom, mass(particle), atol = 10 * ATOL, rtol = RTOL)
    end
end
