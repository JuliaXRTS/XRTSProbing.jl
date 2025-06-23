using QEDprobing
using QEDcore
using QEDprocesses
using Unitful
using BenchmarkTools
using Random
using StatsPlots
using LaTeXStrings


# setup parameter
NE = 1.0e21u"cm^(-3)"
TEMP = 10.0u"keV"

# differential cross section
PROC = Thomson()
MODEL = PerturbativeQED()
PSL = PhotonSphericalLayout(PhotonElectronHeadsOnSystem(XAxis()))
DCS = DifferentialCrossSection(PROC, MODEL, PSL)

# field
FIELD = DistributionBasedField(QEDprobing.GaussianPhotonDist(9.0 * 1000 / QEDprobing.ELECTRONMASS, 0.5 * 1000 / QEDprobing.ELECTRONMASS, XAxis()))

# medium -> non-degenerated approximation
MEDIUM = IdealElectronSystem(NE, TEMP, NoApprox())

# electron dist -> important: same electron temperature
ELEC_DIST = EnergyBasedElectronDistribution(MaxellElectronEnergyDistribution(TEMP))

# kin cuts
CUTS = KinematicCuts(
    (
        1.0,
        7.5 * 1000 / QEDprobing.ELECTRONMASS,
        -1.0,
        0.0,
    ),
    (
        1.2,
        10.5 * 1000 / QEDprobing.ELECTRONMASS,
        1.0,
        2 * pi,
    )
)


# probing setup
probing_setup = ElectronProbingSetup(DCS, FIELD, ELEC_DIST, MEDIUM, CUTS)

# proposal generator
VP = VegasProposal(probing_setup)

## training
RNG = Xoshiro(137)
MAXITER = 100
NCALLS = 10000

train!(RNG, VP, MAXITER, NCALLS)

## max finder
@info "Finding maximum ... "
MAXFINDER = QuantileReductionMethod(0.001, Int(1.0e6))
test_max_weight = QEDprobing.findmax(RNG, probing_setup, MAXFINDER, VP)
@info "Found at: $test_max_weight"

## sampling
N = Int(1.0e7)
BATCH_SIZE = Int(1.0e6)

EG = EventGenerator(probing_setup, VP, test_max_weight)
