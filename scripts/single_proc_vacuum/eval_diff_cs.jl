using QEDprobing
using QEDcore
using QEDprocesses

proc = Thomson() # default: all spin and all pol
mod = PerturbativeQED() # default: vacuum

in_psl = TwoBodyRestSystem() # default: energy of the photon
out_psl = PhotonSphericalLayout(in_psl) # in_coords,cth,phi
# default:     out_psl = PhotonSphericalLayout(in_psl,FixedOmegaPrime()) # omega, cth, phi
# alternative: out_psl = PhotonSphericalLayout(in_psl,FixedOmega()) # omega_prime, cth, phi


# fully open dcs

full_dcs = DifferentialCrossSection(proc, mod, out_psl)

omega = 10.0
cth = range(-1, 1, 100)
phi = 0.1

dcs_vals = map(x -> full_dcs((omega, x, phi)), cth)

# fixed coord dcs
# coord_names(out_psl) # returns: in(:omega), out(:cth,:phi)

fix_coord_dcs = DifferentialCrossSection(proc, mod, out_psl, omega = omega, phi = phi)
dcs_vals = fix_coord_dcs.(cth)
