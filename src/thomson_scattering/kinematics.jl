# TODO: implement PSL
# - for given TwoBodyInPSL the out-PSL just boosts the photon momentum into the elab
# system, builds k_prime and p_prime. And, finally boosts them back into the original
# frame.
# - be aware: if the init electron is in elab already, no boost is necessary!


# Phase Space Definitions

# TODO: rename kinematic mode -> dispersion relation
# kinematic modes
abstract type AbstractKinematicMode end

"""

    is_elatic(::AbstractKinematicMode)

Convenience function to check if a given `AbstractKinematicMode` describes elastic scattering.
"""
function is_elatic end

"""

    Elastic()

Configuration type to ensure, that a given phase space layout describes elatic scattering.
This means, the energy of the beam does not change during scattering.
"""
struct Elastic <: AbstractKinematicMode end
@inline is_elastic(::Elastic) = true

"""

    InElastic()

Configuration type to ensure, that a given phase space layout describes in-elatic scattering.
This means, the energy of the beam does change accordingly during scattering.
"""
struct InElastic <: AbstractKinematicMode end
@inline is_elastic(::InElastic) = false

"""

    ElabPhotonSphSystem()

Represents the ps system with init electron at rest and the out photon described in spherical coordiantes, i.e. polar and azimuthal angle.

"""
struct PhotonSphericalLayout{INPSL <: TwoBodyTargetSystem, K <: AbstractKinematicMode} <:
    AbstractOutPhaseSpaceLayout{INPSL}
    in_psl::INPSL
    kin_mode::K
end

# accessor
kinematic_mode(psl::PhotonSphericalLayout) = psl.kin_mode
is_elastic(psl::PhotonSphericalLayout) = is_elastic(kinematic_mode(psl))

# default
PhotonSphericalLayout(in_psl::TwoBodyTargetSystem) =
    PhotonSphericalLayout(in_psl, InElastic())

# interface
QEDbase.in_phase_space_layout(psl::PhotonSphericalLayout) = psl.in_psl
QEDbase.phase_space_dimension(
    ::AbstractProcessDefinition,
    ::AbstractModelDefinition,
    ::PhotonSphericalLayout,
) = 2

function _TS_omega_prime_elab(::InElastic, om, cth)
    return om / (1 + om * (1 - cth))
end

function _TS_omega_prime_elab(::Elastic, om, cth)
    return om
end

function QEDbase._build_momenta(
        proc::Thomson,
        model::PerturbativeQED,
        in_moms::NTuple{2, <:AbstractFourMomentum},
        psl::PhotonSphericalLayout,
        out_coords::NTuple{2, T},
    ) where {T <: Real}

    P = @inbounds in_moms[1]
    K = @inbounds in_moms[2]
    cth = @inbounds out_coords[1]
    phi = @inbounds out_coords[2]
    om = getE(K)
    omp = _TS_omega_prime_elab(kinematic_mode(psl), om, cth)
    sth = sqrt(1 - cth^2)
    sphi, cphi = sincos(phi)
    Kp = SFourMomentum(omp, omp * sth * cphi, omp * sth * sphi, omp * cth)

    Pp = P + K - Kp
    return Pp, Kp
    moms = _TS_momenta_elab_sph(in_moms, cth, phi)
    return moms
end

function _coordinate_boundaries(::Thomson, ::PerturbativeQED, ::PhotonSphericalLayout)
    return (-1.0, 0.0), (1.0, 2 * pi)
end

function coordinate_symbols(out_psl::PhotonSphericalLayout)
    return (coordinate_symbols(in_phase_space_layout(out_psl))..., :cos_theta, :phi)
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
