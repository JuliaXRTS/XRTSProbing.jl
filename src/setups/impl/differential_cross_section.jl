abstract type AbstractDifferentialCrossSection{P, M, PSL} <: AbstractProcessSetup{P, M, PSL} end

struct DifferentialCrossSection{
        P <: AbstractProcessDefinition,
        M <: AbstractModelDefinition,
        PSL <: AbstractOutPhaseSpaceLayout,
        CC,
        L,
    } <: AbstractDifferentialCrossSection{P, M, PSL}
    proc::P
    model::M
    psl::PSL
    cached_coords::CC

    function DifferentialCrossSection(
            proc::P,
            model::M,
            psl::PSL;
            cached_coords...,
        ) where {
            P <: AbstractProcessDefinition,
            M <: AbstractModelDefinition,
            PSL <: AbstractOutPhaseSpaceLayout,
        }
        L = length(cached_coords)
        # TODO: check if keys(cached_coords) are in coord_syms
        return new{P, M, PSL, typeof(cached_coords), L}(proc, model, psl, cached_coords)
    end
end

# process setup interface

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

coordinate_boundaries(dcs::DifferentialCrossSection) = _coordinate_boundaries(process(dcs), model(dcs), phase_space_layout(dcs))

## Computation

# TODO:
# - this adds an overhead of ~40ns -> needs some optimization!
@inline function _build_dcs_input(diffCS::DifferentialCrossSection, remain_coords::Tuple)
    ccoords = cached_coords(diffCS)
    ccoord_syms = keys(ccoords)
    coord_syms = coordinate_symbols(phase_space_layout(diffCS))
    input = insert_by_name(coord_syms, remain_coords; ccoords...)
    return input
end

# overwrite the generic implementation, to include the cached coords
@inline function _build_psp_with_cached(diffCS::DifferentialCrossSection, remain_coords::Tuple)
    input = _build_dcs_input(diffCS, remain_coords)
    psp = PhaseSpacePoint(
        process(diffCS),
        model(diffCS),
        phase_space_layout(diffCS),
        input
    )

    return psp
end

@inline function _build_psp_without_cached(diffCS::DifferentialCrossSection, coords::Tuple)
    psp = PhaseSpacePoint(
        process(diffCS),
        model(diffCS),
        phase_space_layout(diffCS),
        coords
    )
    return psp
end

# overwrite the generic implementation, to include the cached coords
# everything downstream should work with that, especially `_compute`
function _build_psp(
        diffCS::DifferentialCrossSection{P, M, PS, CC, N},
        coords::Tuple,
    ) where {P, M, PS, CC, N}
    if N == 0
        return _build_psp_without_cached(diffCS, coords)
    else
        return _build_psp_with_cached(diffCS, coords)
    end
end

@inline function _compute(
        diffCS::DifferentialCrossSection{P, M, PS},
        psp::AbstractPhaseSpacePoint{P, M, PS},
    ) where {P, M, PS}
    return differential_cross_section(psp)
end
