using QEDprobing
using QEDcore
using QEDprocesses

using Random
RNG = MersenneTwister(137)

proc = Thomson() # default: all spin and all pol
mod = PerturbativeQED() # default: vacuum

in_psl = TwoBodyRestSystem() # default: energy of the photon
out_psl = PhotonSphericalLayout(in_psl) # in_coord; cth, phi (default: FixedOmegaPrime())

omega = 10.0
phi = 0.1

fix_coord_dcs = DifferentialCrossSection(proc, mod, out_psl, omega = omega, phi = phi)

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
