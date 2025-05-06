using QEDprobing
using QEDcore
using QEDprocesses

proc = Thomson() # default: all spin and all pol

T = 10u"eV"
n_e = 1.0e21u"cm^3"

dsf = LindhardDSF(n_e, T)
mod = PerturbativeQED(dsf) # pertQED within medium described by dsf

# gaussian distributed energy of the second particle with mean 10.0 and std 1.0
in_psl = TwoBodyRestSystem(Gaussian(Energy(2), 10.0, 1.0))

# omega is a function of omega_prime
out_psl = PhotonSphericalLayout(in_psl, FixedOmega()) # omega_prime, cth, phi

phi = 0.1

fix_coord_dcs = DifferentialCrossSection(proc, mod, out_psl, phi = phi) # omega_prime,cth are open, phi is fixed

# vegas proposal
alpha = 0.5

vp = VegasProposal(fix_coord_dcs, alpha)
niter = 100
neval = 1000

train!(RNG, vp, niter, neval)

# max finding
max_finder = QuantileReduction(0.001)
max_weight = findmax(RNG, dcs, max_finder, vp)

# sampling
sampler = EventGenerator(fix_coord_dcs, vp, max_weight)
N = 1.0e6

events = generate_events(RNG, sampler, N)

save_events(events, "test_events.h5")
