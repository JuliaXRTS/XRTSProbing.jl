###############
# patches for QEDprocesses.jl
#
# These functions should be removed if they are in the main line of QEDprocesses.jl
###############

# Opting out the squared matrix elements summed over spins and polarization
# 
# This could also be done for all processes, but here just for (general) Compton to avoid ambiguities
function _matrix_element_square_summed(
    proc::Compton,
    model::AbstractModelDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return sum(
        QEDprocesses._matrix_element_square(proc, model, in_phase_space, out_phase_space)
    )
end

function QEDprocesses._unsafe_differential_probability(
    proc::Compton,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    matrix_elements_sq_summed = _matrix_element_square_summed(
        proc, model, in_phase_space, out_phase_space
    )

    normalization = QEDprocesses._averaging_norm(proc)

    ps_fac = QEDprocesses._phase_space_factor(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )

    return normalization * matrix_elements_sq_summed * ps_fac
end

###
# Special implementation for unpolarized perturbative Compton
#
# Taken from: Peskin, Schroeder "An introduction to Quantum Field Theory" - Eq. (5.105)
###
function _compton_matrix_element_square_spinsum_polsum(in_ps, out_ps)
    pre_fac = 8 * ELEMENTARY_CHARGE^4
    pk = prod(in_ps)
    pkp = in_ps[1] * out_ps[2]

    pkp_over_pk = pkp / pk
    term = 1.0 / (pk) - 1.0 / pkp

    return pre_fac * (pkp_over_pk + 1.0 / pkp_over_pk + 2 * term + term^2)
end

function _matrix_element_square_summed(
    proc::Compton{AllSpin,AllPolarization,AllSpin,AllPolarization},
    model::PerturbativeQED,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _compton_matrix_element_square_spinsum_polsum(in_phase_space, out_phase_space)
end