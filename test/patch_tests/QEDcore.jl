using XRTSProbing
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
INCOORDS = (2.0, 2.0e-2) # E = 2m, omega=10keV
PARTICLE_DIRECTIONS = (XAxis(), YAxis(), ZAxis())


@testset "defaults" begin
    @test PhotonElectronHeadsOnSystem() == PhotonElectronHeadsOnSystem(ZAxis())
end

@testset "dir = $dir" for dir in PARTICLE_DIRECTIONS
    INPSL = PhotonElectronHeadsOnSystem(dir)
    P, K = build_momenta(PROC, MODEL, INPSL, INCOORDS)

    @testset "properties" begin
        @test particle_direction(INPSL) == dir
    end

    @testset "phase space dimension" begin
        @test phase_space_dimension(PROC, MODEL, INPSL) == 2
    end

    @testset "on-shell check" begin
        @test isapprox(getMass(P), mass(Electron()))
        @test isapprox(getMass(K), mass(Photon()))
    end
end
