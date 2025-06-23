using QEDprobing
using QEDcore
using QEDprocesses
using Unitful
using BenchmarkTools
using Random
using Plots


# setup parameter
NE = 1.0e21u"cm^(-3)"
TEMP = 10.0u"keV"
#NE = 2.01e29u"m^(-3)"
#TEMP = 12.53u"eV"


# differential cross section
PROC = Thomson()
MODEL = PerturbativeQED()
PSL = PhotonSphericalLayout(PhotonElectronHeadsOnSystem())
DCS = DifferentialCrossSection(PROC, MODEL, PSL)

# field
FIELD = DistributionBasedField(QEDprobing.GaussianPhotonDist(1.0e-1, 1.0e-3, ZAxis()))
#FIELD = DistributionBasedField(QEDprobing.UniformPhotonDist(2e-3,ZAxis() ))

# medium -> non-degenerated approximation
MEDIUM = IdealElectronSystem(NE, TEMP, NoApprox())

# electron dist -> important: same electron temperature
ELEC_DIST = EnergyBasedElectronDistribution(MaxellElectronEnergyDistribution(TEMP))
@show energy_mean(ELEC_DIST)
@show energy_spectrum(ELEC_DIST, 1.1)
P = plot(
    xscale = :log10,
    xlab = "electron energy [m_e]",
    ylab = "spectrum [a.u.]"
)
plot!(P, x -> energy_spectrum(ELEC_DIST, x), 1.0, 1.2, lab = "energy spectrum")
vline!(P, [energy_mean(ELEC_DIST)], lab = "energy mean")

savefig(P, "elec_dist.pdf")

# kin cuts
CUTS = KinematicCuts(
    (
        1.0,
        1.0e-5,
        -1.0,
        0.0,
    ),
    (
        1.2,
        1.0e-1 + 5 * 1.0e-3,
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
N = Int(1.0e6)
BATCH_SIZE = Int(1.0e5)

EG = EventGenerator(probing_setup, VP, test_max_weight)

@info "generating $N events with batchsize $BATCH_SIZE ..."
event_list = generate_events(RNG, EG, N, BATCH_SIZE)
@info "done."

filename = "events.h5"
QEDprobing.save(event_list, filename)
@info "Saved event list in $filename"

QEDprobing.save_photon_en_mag(event_list, keV = true)
@info "Saved photon energy and magnitude in $filename"

#=
# plotting
momentum_transfer(ev::Event) = momentum_transfer(ev.psp)

function momentum_transfer(psp::AbstractPhaseSpacePoint)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    out_photon_mom = momentum(psp, Outgoing(), Photon())
    return out_photon_mom - in_photon_mom
end

energy_transfer(ev::Event) = energy_transfer(ev.psp)
function energy_transfer(psp::AbstractPhaseSpacePoint)
    return getE(momentum_transfer(psp))
end

three_momentum_transfer(ev::Event) = three_momentum_transfer(ev.psp)
function three_momentum_transfer(psp::AbstractPhaseSpacePoint)
    return getMag(momentum_transfer(psp))
end

in_photon_energy(ev::Event) = in_photon_energy(ev.psp)
function in_photon_energy(psp::AbstractPhaseSpacePoint)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    return getE(in_photon_mom)
end

out_photon_cth(ev::Event) = out_photon_cth(ev.psp)
function out_photon_cth(psp::AbstractPhaseSpacePoint)
    out_photon_mom = momentum(psp, Outgoing(), Photon())
    return getCosTheta(out_photon_mom)
end

om_sample = energy_transfer.(event_list)
q_sample = three_momentum_transfer.(event_list)
omx_sample = in_photon_energy.(event_list)
cth_sample = out_photon_cth.(event_list)

#P = histogram(q_sample,xlim = (0.0,0.01))
P = histogram(q_sample)
savefig(P, "test_q2.pdf")

#P = histogram(om_sample,xlim = (-1e-5,0.0))
P = histogram(om_sample)
savefig(P, "test_om2.pdf")

P = histogram(cth_sample)
savefig(P, "test_cth2.pdf")

P = histogram2d(om_sample, q_sample)
savefig(P, "test2.pdf")

P = histogram(omx_sample)
savefig(P, "test_omx2.pdf")

=#
