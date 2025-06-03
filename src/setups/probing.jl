struct ProbingSetup{T, DOF, D, F, M, C} <: AbstractProcessSetup
    dcs::D
    field::F
    medium::M
    kin_cuts::C

    function ProbingSetup(
            dcs::D,
            field::F,
            medium::M,
            cuts::C,
        ) where {
            N,
            T <: Real,
            D <: DifferentialCrossSection,
            F <: QEDprobing.AbstractSpectrumBasedField,
            M <: AbstractMatterModel, # consider updating to abstract medium?!
            C <: AbstractKinematicCuts{N, T},
        }

        degree_of_freedom(cuts) == degree_of_freedom(dcs) || throw(
            ArgumentError(
                "number of cuts must be equal to the degree-of-freedom given by the differential cross section",
            ),
        )

        return new{T, N, D, F, M, C}(dcs, field, medium, cuts)
    end
end

# accessors
@inline differential_cross_section_setup(stp::ProbingSetup) = stp.dcs
@inline background_field(stp::ProbingSetup) = stp.field
@inline medium(stp::ProbingSetup) = stp.medium
@inline kinematic_cuts(stp::ProbingSetup) = stp.kin_cuts

# interface and deligations
@inline QEDbase.process(stp::ProbingSetup) = process(differential_cross_section_setup(stp))
@inline QEDbase.model(stp::ProbingSetup) = model(differential_cross_section_setup(stp))
@inline QEDbase.phase_space_layout(stp::ProbingSetup) =
    phase_space_layout(differential_cross_section_setup(stp))
@inline QEDbase.in_phase_space_layout(stp::ProbingSetup) =
    in_phase_space_layout(phase_space_layout(stp))
@inline QEDprobing.cached_coords(stp::ProbingSetup) =
    cached_coords(differential_cross_section_setup(stp))

degree_of_freedom(stp::ProbingSetup{T, N}) where {T, N} = N

coordinate_boundaries(stp::ProbingSetup) = extrema(kinematic_cuts(stp))


# compute
function _compute(stp::ProbingSetup, psp::PhaseSpacePoint)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    out_photon_mom = momentum(psp, Outgoing(), Photon())

    #TODO: check order
    mom_transfer = out_photon_mom - in_photon_mom

    in_dist = energy_spectrum(background_field(stp), getE(in_photon_mom))
    #in_dist = 1.0
    # TODO: move prefac to correct position
    dsf = dynamic_structure_factor(medium(stp), (getE(mom_transfer), getMag(mom_transfer)))
    #dsf = 1.0
    hard_dcs = differential_cross_section_setup(stp)(psp)
    #hard_dcs = 1.0

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
    return in_dist * (dsf * hard_dcs)
end

# move this to the generic implementations
function _compute(stp::ProbingSetup{T, N}, coords::NTuple{N, T}) where {T <: Real, N}

    psp = PhaseSpacePoint(process(stp), model(stp), phase_space_layout(stp), coords)
    return _compute(stp, psp)
end

# overwrite the generic one to include check for cuts
function compute(stp::ProbingSetup{T, N}, coords::NTuple{N, T}) where {T <: Real, N}
    all_within_cuts(kinematic_cuts(stp), coords) || return zero(T)

    return _compute(stp, coords)
end

# TODO:
# - rethink and refac this for the sampler!
# - maybe move this to the abstract process setup
@inline function _build_event(
        stp::ProbingSetup,
        coords::Tuple,
        jac::Real,
    )
    psp =
        PhaseSpacePoint(diffCS.proc, diffCS.model, diffCS.psl, coords)
    return Event(psp, _compute(stp, psp) * jac)
end


# TODO:
# - think about implementing a probing weight type
# - this could hold (mom_transfer, in_dist, dsf, hard_dcs) as a result of `compute`
