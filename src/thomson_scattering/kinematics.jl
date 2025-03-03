# TODO: implement PSL
# - for given TwoBodyInPSL the out-PSL just boosts the photon momentum into the elab
# system, builds k_prime and p_prime. And, finally boosts them back into the original
# frame.
# - be aware: if the init electron is in elab already, no boost is necessary!

function _TS_omega_prime_elab(om, cth)
    return om / (1 + om * (1 - cth))
end

function _TS_momenta_elab_sph(in_moms, cth, phi)
    P = @inbounds in_moms[1]
    K = @inbounds in_moms[2]
    om = getE(K)
    omp = _TS_omega_prime_elab(om, cth)
    sth = sqrt(1 - cth^2)
    sphi, cphi = sincos(phi)
    Kp = SFourMomentum(omp, omp * sth * cphi, omp * sth * sphi, omp * cth)

    Pp = P + K - Kp
    return Pp, Kp
end

# Phase Space Definitions

"""

    ElabPhotonSphSystem()

Represents the ps system with init electron at rest and the out photon described in spherical coordiantes, i.e. polar and azimuthal angle.

"""
struct PhotonSphericalLayout{INPSL<:TwoBodyTargetSystem} <:
       AbstractOutPhaseSpaceLayout{INPSL}
    in_psl::INPSL
end

QEDbase.in_phase_space_layout(psl::PhotonSphericalLayout) = psl.in_psl
QEDbase.phase_space_dimension(
    ::AbstractProcessDefinition,
    ::AbstractModelDefinition,
    ::PhotonSphericalLayout,
) = 2


function QEDbase._build_momenta(
    proc::Thomson,
    model::PerturbativeQED,
    in_moms::NTuple{2,<:AbstractFourMomentum},
    ps_def::PhotonSphericalLayout,
    out_coords::NTuple{2,T},
) where {T<:Real}

    cth = @inbounds out_coords[1]
    phi = @inbounds out_coords[2]
    moms = _TS_momenta_elab_sph(in_moms, cth, phi)
    return moms
end

function _coordinate_boundaries(::Thomson, ::PerturbativeQED, ::PhotonSphericalLayout)
    return (-1.0, 0.0), (1.0, 2 * pi)
end

# TODO: extent PSL interface
# - include a function _coordinate_boundaries(proc,model,psl) which return two
# tuples: left and right boundaries of every coordinate.
# - this could also include cuts in the future

# TODO: extent PSL interface
# - degree of freedom == number of coordinates which need to be passed in
# example: CoordMap(Compton,PertQED,PSL) might get (om,cth,phi) -> dof=3
# but: CoordMapCached(Compton,PertQED,PSL,om) gets (cth,phi) -> dof = 2

# For now, we work on PSDEFs, cached DiffCS and simple implementations of the
# boundaries
