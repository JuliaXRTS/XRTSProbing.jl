abstract type AbstractProbingSetup{P, M, PSL} <: AbstractProcessSetup{P, M, PSL} end

struct ElectronProbingSetup{T, DOF, P, M, PSL, D, F, E, MEDIUM, C} <: AbstractProbingSetup{P, M, PSL}
    dcs::D
    field::F
    elec_dist::E
    medium::MEDIUM
    kin_cuts::C

    function ElectronProbingSetup(
            dcs::D,
            field::F,
            elec_dist::E,
            medium::MEDIUM,
            cuts::C,
        ) where {
            N,
            T <: Real,
            P <: AbstractProcessDefinition,
            M <: AbstractModelDefinition,
            PSL <: AbstractOutPhaseSpaceLayout,
            D <: AbstractDifferentialCrossSection{P, M, PSL},
            F <: AbstractSpectrumBasedField,
            E <: AbstractElectronDistribution,
            MEDIUM <: AbstractMatterModel, # consider updating to abstract medium?!
            C <: AbstractKinematicCuts{N, T},
        }

        # TODO: remove this and introduce (N,T) as type parameter to AbstractProcessSetup
        # (and therefore AbstractDifferentialCrossSection)
        degree_of_freedom(cuts) == degree_of_freedom(dcs) || throw(
            ArgumentError(
                "number of cuts must be equal to the degree-of-freedom given by the differential cross section",
            ),
        )

        return new{T, N, P, M, PSL, D, F, E, MEDIUM, C}(dcs, field, elec_dist, medium, cuts)
    end
end

# accessors
@inline differential_cross_section_setup(stp::ElectronProbingSetup) = stp.dcs
@inline background_field(stp::ElectronProbingSetup) = stp.field
@inline electron_distribution(stp::ElectronProbingSetup) = stp.elec_dist
@inline medium(stp::ElectronProbingSetup) = stp.medium
@inline kinematic_cuts(stp::ElectronProbingSetup) = stp.kin_cuts

# interface and deligations
@inline QEDbase.process(stp::ElectronProbingSetup) = process(differential_cross_section_setup(stp))
@inline QEDbase.model(stp::ElectronProbingSetup) = model(differential_cross_section_setup(stp))
@inline QEDbase.phase_space_layout(stp::ElectronProbingSetup) =
    phase_space_layout(differential_cross_section_setup(stp))
@inline QEDbase.in_phase_space_layout(stp::ElectronProbingSetup) =
    in_phase_space_layout(phase_space_layout(stp))
@inline QEDprobing.cached_coords(stp::ElectronProbingSetup) =
    cached_coords(differential_cross_section_setup(stp))

degree_of_freedom(stp::ElectronProbingSetup{T, N}) where {T, N} = N

coordinate_boundaries(stp::ElectronProbingSetup) = extrema(kinematic_cuts(stp))

# compute
function _compute(stp::ElectronProbingSetup, psp::PhaseSpacePoint)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    in_elec_mom = momentum(psp, Incoming(), Electron())
    out_photon_mom = momentum(psp, Outgoing(), Photon())

    #TODO: check order
    mom_transfer = out_photon_mom - in_photon_mom

    in_dist = energy_spectrum(background_field(stp), getE(in_photon_mom))

    elec_dist = energy_spectrum(electron_distribution(stp), getE(in_elec_mom))

    # TODO: move prefac to correct position
    dsf = dynamic_structure_factor(medium(stp), (getE(mom_transfer), getMag(mom_transfer)))

    hard_dcs = differential_cross_section_setup(stp)(psp)

    #=
    if isinf(dsf) || isinf(in_dist) || isinf(hard_dcs) || isnan(dsf) || isnan(in_dist) || isnan(hard_dcs)
        throw(
            InvalidInputError(
                """
                something went wrong at om = $(getE(mom_transfer)) q = $(getMag(mom_transfer))
                with
                dsf = $dsf
                in_dist= $in_dist
                hard_dcs = $hard_dcs
                """
            )
        )
    end
    =#
    return elec_dist * (in_dist * (dsf * hard_dcs))
end

# overwrite generic fallback to include dof check in the types
# TODO: remove this, after it is brought upstream
function _compute(stp::ElectronProbingSetup{T, N}, coords::NTuple{N, T}) where {T <: Real, N}
    psp = PhaseSpacePoint(process(stp), model(stp), phase_space_layout(stp), coords)
    return _compute(stp, psp)
end

# overwrite the generic one to include check for cuts
function compute(stp::ElectronProbingSetup{T, N}, coords::NTuple{N, T}) where {T <: Real, N}
    all_within_cuts(kinematic_cuts(stp), coords) || return zero(T)

    return _compute(stp, coords)
end
