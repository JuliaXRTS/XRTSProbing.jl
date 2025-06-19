using QEDprobing
using QEDcore
using QEDprocesses
using Random

ATOL = 1.0e-12
RTOL = 1.0e-6

RNG = MersenneTwister(137137)

IN_PARTICLES = (Electron(), Photon())
OUT_PARTICLES = (Electron(), Photon())

PROC = ScatteringProcess(IN_PARTICLES, OUT_PARTICLES)

MODEL = PertQED()
INPSL = PhotonElectronHeadsOnSystem()
INCOORDS = (2.0, 2.0e-2) # E = 2m, omega=10keV

P, K = build_momenta(PROC, MODEL, INPSL, INCOORDS)

@testset "phase space dimension" begin
    @test phase_space_dimension(PROC, MODEL, INPSL) == 2
end

@testset "on-shell check" begin
    @test isapprox(getMass(P), mass(Electron()))
    @test isapprox(getMass(K), mass(Photon()))
end
