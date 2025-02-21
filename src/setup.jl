
abstract type AbstractDifferentialCrossSection end

struct DifferentialCrossSection{
    PROC<:AbstractProcessDefinition,
    MODEL<:AbstractModelDefinition,
    PSDEF<:AbstractPhasespaceDefinition,
} <: AbstractDifferentialCrossSection
    proc::PROC
    model::MODEL
    ps_def::PSDEF
end

QEDcore.process(d::DifferentialCrossSection) = d.proc
QEDcore.model(d::DifferentialCrossSection) = d.model
QEDcore.phase_space_definition(d::DifferentialCrossSection) = d.ps_def

function (diffCS::DifferentialCrossSection)(in_coords::Tuple, out_coords::Tuple)
    psp = PhaseSpacePoint(diffCS.proc, diffCS.model, diffCS.ps_def, in_coords, out_coords)
    return differential_cross_section(psp)
end

function (diffCS::DifferentialCrossSection{P,M,PS})(
    psp::AbstractPhaseSpacePoint{P,M,PS},
) where {P,M,PS}
    return differential_cross_section(psp)
end

struct DifferentialCrossSectionCached{
    PROC<:AbstractProcessDefinition,
    MODEL<:AbstractModelDefinition,
    PSDEF<:AbstractPhasespaceDefinition,
    COORDS<:Tuple,
} <: AbstractDifferentialCrossSection
    proc::PROC
    model::MODEL
    ps_def::PSDEF
    in_coords::COORDS
end
Base.broadcastable(d::DifferentialCrossSectionCached) = Ref(d)

QEDcore.process(d::DifferentialCrossSectionCached) = d.proc
QEDcore.model(d::DifferentialCrossSectionCached) = d.model
QEDcore.phase_space_definition(d::DifferentialCrossSectionCached) = d.ps_def
in_coordinates(d::DifferentialCrossSectionCached) = d.in_coords

function (diffCS::DifferentialCrossSectionCached)(out_coords::Tuple)
    psp = PhaseSpacePoint(
        diffCS.proc,
        diffCS.model,
        diffCS.ps_def,
        diffCS.in_coords,
        out_coords,
    )
    return differential_cross_section(psp)
end
function (diffCS::DifferentialCrossSectionCached{P,M,PS})(
    psp::AbstractPhaseSpacePoint{P,M,PS},
) where {P,M,PS}
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
    psp = PhaseSpacePoint(
        diffCS.proc,
        diffCS.model,
        diffCS.ps_def,
        diffCS.in_coords,
        out_coords,
    )
    return Event(psp, differential_cross_section(psp) * jac)
end

@inline function _build_event(diffCS::DifferentialCrossSectionCached, out_coords::Tuple)
    psp = PhaseSpacePoint(
        diffCS.proc,
        diffCS.model,
        diffCS.ps_def,
        diffCS.in_coords,
        out_coords,
    )
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
