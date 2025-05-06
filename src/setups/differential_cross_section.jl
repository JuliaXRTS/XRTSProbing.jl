abstract type AbstractDifferentialCrossSection <: AbstractProcessSetup end

struct DifferentialCrossSection{
        PROC <: AbstractProcessDefinition,
        MODEL <: AbstractModelDefinition,
        PSL <: AbstractOutPhaseSpaceLayout,
        CC,
        L,
    } <: AbstractDifferentialCrossSection
    proc::PROC
    model::MODEL
    psl::PSL
    cached_coords::CC

    function DifferentialCrossSection(
            proc::PROC,
            model::MODEL,
            psl::PSL;
            cached_coords...,
        ) where {
            PROC <: AbstractProcessDefinition,
            MODEL <: AbstractModelDefinition,
            PSL <: AbstractOutPhaseSpaceLayout,
        }
        L = length(cached_coords)
        # TODO: check if keys(cached_coords) are in coord_syms
        return new{PROC, MODEL, PSL, typeof(cached_coords), L}(proc, model, psl, cached_coords)
    end
end

@inline QEDbase.process(d::DifferentialCrossSection) = d.proc
@inline QEDbase.model(d::DifferentialCrossSection) = d.model
@inline QEDbase.phase_space_layout(d::DifferentialCrossSection) = d.psl
@inline QEDbase.in_phase_space_layout(d::DifferentialCrossSection) =
    in_phase_space_layout(d.psl)
@inline cached_coords(d::DifferentialCrossSection) = d.cached_coords
@inline function degree_of_freedom(dcs::DifferentialCrossSection)
    in_dim = phase_space_dimension(process(dcs), model(dcs), in_phase_space_layout(dcs))
    out_dim = phase_space_dimension(process(dcs), model(dcs), phase_space_layout(dcs))
    n_cached = length(cached_coords(dcs))
    return in_dim + out_dim - n_cached
end

# TODO: remove that
@inline function _assert_valid_input(dcs::DifferentialCrossSection, input::Tuple)
    n_input = length(input)
    degree_of_freedom(dcs) == n_input || throw(
        DimensionMismatch(
            """
            number of input coordinates must match degree of freedom of the differential cross section

            expected: degree_of_freedom(dcs) == $(degree_of_freedom(dcs))
            given: length(input) == $(n_input)
            """,
        ),
    )
    return nothing
end

# TODO:
# - this adds an overhead of ~40ns -> needs some optimization!
@inline function _compute_with_cached(
        diffCS::DifferentialCrossSection,
        remain_coords::Tuple,
    )
    ccoords = cached_coords(diffCS)
    ccoord_syms = keys(ccoords)
    coord_syms = coordinate_symbols(phase_space_layout(diffCS))
    input = insert_by_name(coord_syms, remain_coords; ccoords...)

    return _compute_without_cached(diffCS, input)
end

# TODO: input validation?
@inline function _compute_without_cached(diffCS::DifferentialCrossSection, coords::Tuple)
    in_dim =
        phase_space_dimension(process(diffCS), model(diffCS), in_phase_space_layout(diffCS))
    out_dim =
        phase_space_dimension(process(diffCS), model(diffCS), phase_space_layout(diffCS))
    return _compute(
        diffCS,
        ntuple(i -> coords[i], in_dim),
        ntuple(i -> coords[in_dim + i], out_dim),
    )
end

function _compute_without_cached_barrier(diffCS, coords)
    return _compute_without_cached(diffCS, coords)
end

@inline function _compute(
        diffCS::DifferentialCrossSection{P, M, PS, CC, N},
        coords::Tuple,
    ) where {P, M, PS, CC, N}
    if N == 0
        return _compute_without_cached(diffCS, coords)
    else
        return _compute_with_cached(diffCS, coords)
    end
end

@inline function _compute(
        diffCS::DifferentialCrossSection,
        in_coords::Tuple,
        out_coords::Tuple,
    )
    psp = PhaseSpacePoint(diffCS.proc, diffCS.model, diffCS.psl, in_coords, out_coords)
    return differential_cross_section(psp)
end

@inline function _compute(
        diffCS::DifferentialCrossSection{P, M, PS},
        psp::AbstractPhaseSpacePoint{P, M, PS},
    ) where {P, M, PS}
    return differential_cross_section(psp)
end

# TODO: generalize this call to AbstractSetup
@inline function (diffCS::DifferentialCrossSection{P, M, PS, CC, N})(
        coords::Tuple,
    ) where {P, M, PS, CC, N}
    return compute(diffCS, coords)
end

# TODO: generalize this call to AbstractSetup
@inline function (diffCS::DifferentialCrossSection)(in_coords::Tuple, out_coords::Tuple)
    return compute(diffCS, (in_coords..., out_coords...))
end

# TODO: generalize this call to AbstractProcessSetup
@inline function (diffCS::DifferentialCrossSection{P, M, PS})(
        psp::AbstractPhaseSpacePoint{P, M, PS},
    ) where {P, M, PS}
    return differential_cross_section(psp)
end


# TODO: remove that
struct DifferentialCrossSectionCached{
        PROC <: AbstractProcessDefinition,
        MODEL <: AbstractModelDefinition,
        PSL <: AbstractOutPhaseSpaceLayout,
        COORDS <: Tuple,
    } <: AbstractDifferentialCrossSection
    proc::PROC
    model::MODEL
    psl::PSL
    in_coords::COORDS
end
Base.broadcastable(d::DifferentialCrossSectionCached) = Ref(d)

QEDcore.process(d::DifferentialCrossSectionCached) = d.proc
QEDcore.model(d::DifferentialCrossSectionCached) = d.model
QEDcore.phase_space_layout(d::DifferentialCrossSectionCached) = d.psl
in_coordinates(d::DifferentialCrossSectionCached) = d.in_coords

function (diffCS::DifferentialCrossSectionCached)(out_coords::Tuple)
    psp =
        PhaseSpacePoint(diffCS.proc, diffCS.model, diffCS.psl, diffCS.in_coords, out_coords)
    return differential_cross_section(psp)
end
function (diffCS::DifferentialCrossSectionCached{P, M, PS})(
        psp::AbstractPhaseSpacePoint{P, M, PS},
    ) where {P, M, PS}
    # note: the in_coords of the psp might differ from the one stored in diffCS.
    # However this should not be a problem, because either the psp is constructed using
    # diffCS.in_coords or the user just wants the diff. CS for the given psp. Both is fine
    # and we should not complain.
    return differential_cross_section(psp)
end

@inline function _build_event(
        diffCS::DifferentialCrossSectionCached,
        out_coords::Tuple,
        jac::Real,
    )
    psp =
        PhaseSpacePoint(diffCS.proc, diffCS.model, diffCS.psl, diffCS.in_coords, out_coords)
    return Event(psp, differential_cross_section(psp) * jac)
end

@inline function _build_event(diffCS::DifferentialCrossSectionCached, out_coords::Tuple)
    psp =
        PhaseSpacePoint(diffCS.proc, diffCS.model, diffCS.psl, diffCS.in_coords, out_coords)
    return Event(psp, differential_cross_section(psp))
end

# TODO: special implementations
# - if the diff cs is known directly from coordinates, there is no need to construct the
# psp in the first place
# - consider implementing diff_cs setups for these special cases
# - implement a call on psp (ensure matching proc, model and psdef at compile time); this
# allows for broadcasting the dcs call over arrays of psp, e.g. for GPU execution
# - consider inrolling the arguments (e.g. diff_cs(x,y,z) = diff_cs((x,y,z))), maybe makes
# it easier to call, or to map, or to broadcast.

# TODO: full event
# - consider implementing a function producing (psp,diff_cs), i.e. the actual event
# - but keep the function only returning diff_cs (mostly used for training)
# - consider to revise the actual setup interface (with validate, compute, ... )
