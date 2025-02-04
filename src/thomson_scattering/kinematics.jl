# TODO: implement PSL
# - for given TwoBodyInPSL the out-PSL just boosts the photon momentum into the elab
# system, builds k_prime and p_prime. And, finally boosts them back into the original
# frame.
# - be aware: if the init electron is in elab already, no boost is necessary!

function _TS_omega_prime_elab(om, cth)
    return om / (1 + om * (1 - cth))
end

function _TS_init_electron_mom_elab(m)
    return SFourMomentum(m, zero(m), zero(m), zero(m))
end

function _TS_init_photon_mom_zaligned(om)
    return SFourMomentum(om, zero(om), zero(om), om)
end

function _TS_momenta_elab_sph(om, cth, phi)
    P = _TS_init_electron_mom_elab(one(om))
    K = _TS_init_photon_mom_zaligned(om)
    omp = _TS_omega_prime_elab(om, cth)
    sth = sqrt(1 - cth^2)
    sphi, cphi = sincos(phi)
    Kp = SFourMomentum(omp, omp * sth * cphi, omp * sth * sphi, omp * cth)

    Pp = P + K - Kp
    return P, K, Pp, Kp
end

# Phase Space Definitions

"""

    ElabPhotonSphSystem()

Represents the ps system with init electron at rest and the out photon described in spherical coordiantes, i.e. polar and azimuthal angle.

"""
struct ElabPhotonSphSystem <: AbstractPhasespaceDefinition end

function QEDbase._generate_incoming_momenta(
    proc::ThomsonScattering,
    model::PerturbativeQED,
    in_ps_def::ElabPhotonSphSystem,
    in_ps::NTuple{1,T},
) where {T<:Real}
    om = @inbounds in_ps[1]

    P = _TS_init_electron_mom_elab(one(om))
    K = _TS_init_photon_mom_zaligned(om)

    return (P, K)
end

function QEDbase._generate_outgoing_momenta(
    proc::ThomsonScattering,
    model::PerturbativeQED,
    phase_space_def::ElabPhotonSphSystem,
    in_phase_space::NTuple{1,T},
    out_phase_space::NTuple{2,T},
) where {T<:Real}
    om = @inbounds in_phase_space[1]
    cth = @inbounds out_phase_space[1]
    phi = @inbounds out_phase_space[2]
    moms = _TS_momenta_elab_sph(om, cth, phi)
    return @inbounds (moms[3], moms[4])
end
