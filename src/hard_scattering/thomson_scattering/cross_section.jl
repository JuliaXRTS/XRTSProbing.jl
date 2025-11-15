# General Thomson scattering

@inline function QEDbase._incident_flux(
        psp::InPhaseSpacePoint{<:Thomson, PerturbativeQED, <:PhotonSphericalLayout},
    )
    om = getE(momentum(psp, Incoming(), Photon()))
    return om
end

@inline function QEDbase._averaging_norm(::Type{T}, proc::Thomson) where {T <: Real}
    return inv(T(incoming_multiplicity(proc)))
end

# TODO: consider putting on-shell check for inelastic scattering in here
@inline function QEDbase._is_in_phasespace(psp::PhaseSpacePoint{<:Thomson, PerturbativeQED})
    return isapprox(sum(momenta(psp, Incoming())), sum(momenta(psp, Outgoing())))
end

@inline function QEDbase._phase_space_factor(
        psp::PhaseSpacePoint{<:Thomson, PerturbativeQED, <:PhotonSphericalLayout},
    )
    om = getE(momentum(psp, Incoming(), Photon()))
    return om / (16 * pi^2)
end

# polarized
function QEDbase._matrix_element_square(
        psp::PhaseSpacePoint{
            Thomson{AllSpin, INPOL, AllSpin, OUTPOL},
            PerturbativeQED,
            <:PhotonSphericalLayout,
        },
    ) where {INPOL <: AbstractDefinitePolarization, OUTPOL <: AbstractDefinitePolarization}
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
        psp::PhaseSpacePoint{PROC, PerturbativeQED, <:PhotonSphericalLayout},
    ) where {PROC <: Thomson{AllSpin, AllPol, AllSpin, AllPol}}
    cth = getCosTheta(momentum(psp, Outgoing(), Photon()))
    return 8 * ELEMENTARY_CHARGE_SQUARE^2 * (one(cth) + cth^2)
end

# return the cos of the angle between the given momenta
_cos_delta_theta(mom1, mom2) = dot(view(mom1, 2:4), view(mom2, 2:4)) / (getMag(mom1) * getMag(mom2))

function QEDbase.unsafe_differential_cross_section(
        psp::PhaseSpacePoint{
            Thomson{AllSpin, AllPol, AllSpin, AllPol},
            PerturbativeQED,
            <:PhotonSphericalLayout,
        },
    )
    k = momentum(psp, Incoming(), Photon())
    k_prime = momentum(psp, Outgoing(), Photon())
    cth = _cos_delta_theta(k, k_prime)
    return XRTSProbing._TS_diffCS_pol_spin_sum(cth)
end

function _TS_diffCS_pol_spin_sum(cth)
    return ALPHA_SQUARE / 2 * (1 + cth^2)
end
