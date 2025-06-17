const HBARC = 197.3269718 # MeV fm
const HBARC_eV_ANG = HBARC * 1.0e6 * 1.0e-5 # eV Ang
const ELECTRONMASS = 0.5109989461e6 # eV
const ALPHA = 1 / (137.035999074)
const ALPHA_SQUARE = ALPHA^2
const ELEMENTARY_CHARGE_SQUARED = 4 * pi * ALPHA
const ELEMENTARY_CHARGE = sqrt(ELEMENTARY_CHARGE_SQUARED)
const BOHR_RADIUS_ANG = 0.529177210544 # Ang
const HARTREE = 27.211386245981 # eV

# simple unit transforms
_bohr2ang(x_bohr, p = 1) = x_bohr * (BOHR_RADIUS_ANG^p)
_bohr2cm(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG * 1.0e-8)^p)
_bohr2fm(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG * 1.0e5)^p)
_bohr2inv_eV(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG / HBARC_eV_ANG)^p)
_bohr2inv_me(x_bohr, p = 1) = x_bohr * ((BOHR_RADIUS_ANG / HBARC_eV_ANG * ELECTRONMASS)^p)

_inv_bohr2ang(x_inv_bohr, p = 1) = x_inv_bohr / (BOHR_RADIUS_ANG^p)
_inv_bohr2cm(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG * 1.0e-8)^p)
_inv_bohr2fm(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG * 1.0e5)^p)
_inv_bohr2eV(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG / HBARC_eV_ANG)^p)
_inv_bohr2me(x_inv_bohr, p = 1) = x_inv_bohr / ((BOHR_RADIUS_ANG / HBARC_eV_ANG * ELECTRONMASS)^p)

_ang2inv_eV(x_ang, p = 1) = x_ang / (HBARC_eV_ANG^p)
_inv_ang2eV(x_inv_ang, p = 1) = x_inv_ang * (HBARC_eV_ANG^p)

_me2eV(x_me, p = 1) = x_me * (ELECTRONMASS)^p
_me2keV(x_me, p = 1) = x_me * (ELECTRONMASS * 1.0e-3)^p
_me2MeV(x_me, p = 1) = x_me * (ELECTRONMASS * 1.0e-6)^p
_me2hartree(x_me, p = 1) = x_me * (ELECTRONMASS / HARTREE)^p
_me2inv_ang(x_me, p = 1) = x_me * (ELECTRONMASS / HBARC_eV_ANG)^p
_me2inv_bohr(x_me, p = 1) = x_me * (ELECTRONMASS / HBARC_eV_ANG * BOHR_RADIUS_ANG)^p
