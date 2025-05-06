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

# omega_prime is a function of omega
# out_psl = PhotonSphericalLayout(in_psl,FixedOmegaPrime()) # omega, cth, phi

# fully open dcs

full_dcs = DifferentialCrossSection(proc, mod, out_psl)

omega_prime = 8.0
cth = range(-1, 1, 100)
phi = 0.1

dcs_vals = map(x -> full_dcs((omega_prime, x, phi)), cth)

# fixed coord dcs
# coord_names(out_psl) # returns: in(:omega), out(:cth,:phi)

fix_coord_dcs =
    DifferentialCrossSection(proc, mod, out_psl, omega_prime = omega_prime, phi = phi)
dcs_vals = fix_coord_dcs.(cth)
