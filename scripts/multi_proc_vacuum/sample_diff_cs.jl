using QEDprobing
using QEDcore
using QEDprocesses

using Random
RNG = MersenneTwister(137)


# process 1: Thomson
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
sampler1 = EventGenerator(fix_coord_dcs, vp, max_weight)



# process 2: trident
proc = Trident() # default: all spin and all pol
mod = PerturbativeQED() # default: vacuum

in_psl = TwoBodyRestSystem() # default: energy of the photon
out_psl = TridentSphericalSystem(in_psl) # in_coord; cth, phi (default: FixedOmegaPrime())

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
sampler2 = EventGenerator(fix_coord_dcs, vp, max_weight)



## Combined event generation
# two step sampling
# 1. binomial sample based on totCS: choose which process
# 2. standard sampling of the choosen process

combined_sampler = CombinedEventGenerator(sampler1, sampler2)

N = 1e6

events = generate_events(RNG, combined_sampler, N)

save_events(events, "test_events.h5")
