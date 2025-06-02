using QEDprobing
using QEDcore
using QEDprocesses
using Unitful
using BenchmarkTools
using Random
using Plots

# differential cross section
PROC = Thomson()
MODEL = PerturbativeQED()
PSL = PhotonSphericalLayout(TwoBodyTargetSystem())
DCS = DifferentialCrossSection(PROC, MODEL, PSL)

# field
FIELD = DistributionBasedField(QEDprobing.GaussianPhotonDist(1.0e-3, 1.0e-3, ZAxis()))
#FIELD = DistributionBasedField(QEDprobing.UniformPhotonDist(2e-3,ZAxis() ))

# medium
MEDIUM = IdealElectronSystem(2.01e29u"m^(-3)", 12.53u"eV")

# kin cuts
CUTS = KinematicCuts(
    (
        1.0e-6,
        -1.0,
        0.0,
    ),
    (
        2.0e-3,
        1.0,
        2 * pi,
    )
)

# probing setup
probing_setup = ProbingSetup(DCS, FIELD, MEDIUM, CUTS)

# proposal generator
VP = VegasProposal(probing_setup)

## training
RNG = Xoshiro(137)
MAXITER = 100
NCALLS = 10000

train!(RNG, VP, MAXITER, NCALLS)

## sampling
N = Int(1.0e6)
sample, jac = QEDprobing._generate_coords(RNG, VP, N)
sample_psps = QEDprobing._build_psp.(probing_setup, sample)

# plotting
function momentum_transfer(psp)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    out_photon_mom = momentum(psp, Outgoing(), Photon())
    return out_photon_mom - in_photon_mom
end

function energy_transfer(psp)
    return getE(momentum_transfer(psp))
end

function three_momentum_transfer(psp)
    return getMag(momentum_transfer(psp))
end

function in_photon_energy(psp)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    return getE(in_photon_mom)
end

function out_photon_cth(psp)
    out_photon_mom = momentum(psp, Outgoing(), Photon())
    return getCosTheta(out_photon_mom)
end

om_sample = energy_transfer.(sample_psps)
q_sample = three_momentum_transfer.(sample_psps)
omx_sample = in_photon_energy.(sample_psps)
cth_sample = out_photon_cth.(sample_psps)

#P = histogram(q_sample,xlim = (0.0,0.01))
P = histogram(q_sample)
savefig(P, "test_q.pdf")

#P = histogram(om_sample,xlim = (-1e-5,0.0))
P = histogram(om_sample)
savefig(P, "test_om.pdf")

P = histogram(cth_sample)
savefig(P, "test_cth.pdf")

P = histogram2d(om_sample, q_sample)
savefig(P, "test.pdf")

P = histogram(omx_sample)
savefig(P, "test_omx.pdf")
