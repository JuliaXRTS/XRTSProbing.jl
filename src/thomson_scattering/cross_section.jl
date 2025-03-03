
# General Thomson scattering

@inline function QEDbase._incident_flux(
    psp::InPhaseSpacePoint{<:Thomson,PerturbativeQED,<:PhotonSphericalLayout},
)
    om = getE(momentum(psp, Incoming(), Photon()))
    return om
end

@inline function QEDbase._averaging_norm(proc::Thomson)
    return inv(incoming_multiplicity(proc))
end

@inline function QEDbase._is_in_phasespace(psp::PhaseSpacePoint{<:Thomson,PerturbativeQED})
    return isapprox(sum(momenta(psp, Incoming())), sum(momenta(psp, Outgoing())))
end

@inline function QEDbase._phase_space_factor(
    psp::PhaseSpacePoint{<:Thomson,PerturbativeQED,<:PhotonSphericalLayout},
)
    om = getE(momentum(psp, Incoming(), Photon()))
    return om / (16 * pi^2)
end

# polarized
function QEDbase._matrix_element_square(
    psp::PhaseSpacePoint{
        Thomson{AllSpin,INPOL,AllSpin,OUTPOL},
        PerturbativeQED,
        <:PhotonSphericalLayout,
    },
) where {INPOL<:AbstractDefinitePolarization,OUTPOL<:AbstractDefinitePolarization}
    in_photon_mom = momentum(psp, Incoming(), Photon())
    in_photon_pol = _spin_or_pol(process(psp), Photon(), Incoming())
    in_pol_vec = base_state(Photon(), Incoming(), in_photon_mom, in_photon_pol)


    out_photon_mom = momentum(psp, Outgoing(), Photon())
    out_photon_pol = _spin_or_pol(process(psp), Photon(), Outgoing())
    out_pol_vec = base_state(Photon(), Outgoing(), out_photon_mom, out_photon_pol)

    pol_prod = in_pol_vec * out_pol_vec


    return 8 * ELEMENTARY_CHARGE_SQUARE^2 * pol_prod^2
end

# unpolarized
function QEDbase._matrix_element_square(
    psp::PhaseSpacePoint{PROC,PerturbativeQED,<:PhotonSphericalLayout},
) where {PROC<:Thomson{AllSpin,AllPol,AllSpin,AllPol}}
    cth = getCosTheta(momentum(psp, Outgoing(), Photon()))
    return 8 * ELEMENTARY_CHARGE_SQUARE^2 * (one(cth) + cth^2)
end

function QEDbase.unsafe_differential_cross_section(
    psp::PhaseSpacePoint{
        Thomson{AllSpin,AllPol,AllSpin,AllPol},
        PerturbativeQED,
        <:PhotonSphericalLayout,
    },
)
    k_prime = momentum(psp, Outgoing(), Photon())
    cth = getCosTheta(k_prime)
    return QEDprobing._TS_diffCS_pol_spin_sum(cth)
end

function _TS_diffCS_pol_spin_sum(cth)
    return ALPHA_SQUARE / 2 * (1 + cth^2)
end
