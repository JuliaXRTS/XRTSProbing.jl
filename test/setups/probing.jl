using Test
using Random
using Unitful
using Distributions

using QEDcore
using QEDbase
using QEDprocesses
using QEDprobing

RNG = Xoshiro(137)
ATOL = sqrt(eps())
RTOL = sqrt(eps())

# DCS
const PROC = Thomson()
const MODEL = PerturbativeQED()
const INPSL = TwoBodyTargetSystem()
const OUTPSL = PhotonSphericalLayout(INPSL)
const DCS = DifferentialCrossSection(PROC, MODEL, OUTPSL)

# FIELD
const OMX = 2.0e-2 # == 10keV
const OMX_STD = 1.0e-3 # == 0.5 keV
const PHOTONAXIS = ZAxis()
const PHOTONDIST = GaussianPhotonDist(OMX, OMX_STD, PHOTONAXIS)
const FIELD = DistributionBasedField(PHOTONDIST)

# MEDIUM
const NE = 1.0e21u"cm^(-3)"
const MATTERMODEL = IdealElectronSystem{ZeroTemperature}(NE)
const SCREENING = Screening()
const MEDIUM = InteractingElectronSystem(MATTERMODEL, SCREENING)

# CUTS
LOW_OM, UP_OM = 1.0e-2 - 5.0e-3, 1.0e-2 + 5.0e-3
const CUTS = KinematicCuts(
    (
        LOW_OM, # min omx
        -1.0, # min cos theta
        0.0, # min phi
    ), (
        UP_OM, # max omx
        1.0, # max cos theta
        2 * pi, # max phi
    )
)

TEST_PROBING_SETUP = ProbingSetup(DCS, FIELD, MEDIUM, CUTS)

# valid coordinates
const OMEGAS = (rand(RNG, Uniform(LOW_OM, UP_OM)))

# TODO: check cth=1.0 (currently, dsf->NaN)
const CTHS = (-1.0, 0.0, rand(RNG), -rand(RNG))
const PHIS = (0.0, 2 * pi, rand(RNG) * 2 * pi)
const COORDS = (OMEGAS, CTHS, PHIS)


@testset "accessor" begin
    @test differential_cross_section_setup(TEST_PROBING_SETUP) == DCS
    @test process(TEST_PROBING_SETUP) == PROC
    @test model(TEST_PROBING_SETUP) == MODEL
    @test phase_space_layout(TEST_PROBING_SETUP) == OUTPSL
    @test in_phase_space_layout(TEST_PROBING_SETUP) == INPSL

    @test medium(TEST_PROBING_SETUP) == MEDIUM
    @test background_field(TEST_PROBING_SETUP) == FIELD
    @test kinematic_cuts(TEST_PROBING_SETUP) == CUTS
    @test degree_of_freedom(TEST_PROBING_SETUP) == 3
end

@testset "cth = $cth, phi = $phi" for (cth, phi) in Iterators.product(CTHS, PHIS)

    @testset "sanity check" begin
        @testset "om = $om" for om in OMEGAS
            psp = PhaseSpacePoint(PROC, MODEL, OUTPSL, (om,), (cth, phi))
            dcs = DCS(psp)
            in_photon_mom = momentum(psp, Incoming(), Photon())
            out_photon_mom = momentum(psp, Outgoing(), Photon())

            #TODO: check order
            mom_transfer = out_photon_mom - in_photon_mom
            dsf =
                dynamic_structure_factor(MEDIUM, (getE(mom_transfer), getMag(mom_transfer)))
            ph_dist = energy_spectrum(FIELD, om)

            @test isapprox(
                QEDprobing._compute(TEST_PROBING_SETUP, (om, cth, phi)),
                ph_dist * dsf * dcs,
            )

            @test QEDprobing._compute(TEST_PROBING_SETUP, (om + UP_OM, cth, phi)) ==
                zero(om)
        end
    end

end
