###############
# The process setup
#
# In this file, we define the interface for general computation and process setups.
###############

abstract type AbstractComputationSetup end

@inline function _assert_valid_input(stp::AbstractComputationSetup, input)
    return nothing
end

@inline function _post_processing(stp::AbstractComputationSetup, input, result)
    return result
end

function _compute end


# TODO: extent list of derived/delegated accessors
function compute(stp::AbstractComputationSetup, input)
    _assert_valid_input(stp, input)
    raw_result = _compute(stp, input)
    return _post_processing(stp, input, raw_result)
end

abstract type AbstractProcessSetup <: AbstractComputationSetup end


@inline QEDbase.number_incoming_particles(stp::AbstractProcessSetup) =
    number_incoming_particles(process(stp))
@inline QEDbase.number_outgoing_particles(stp::AbstractProcessSetup) =
    number_outgoing_particles(process(stp))

# TODO:
# - add `coordinate_boundaries` to the interface (dcs-> infer from psl, pstp -> cuts)
# - use `compute` to check if coords are in bounds (use _compute if the check does not to
# be done)
# - don't check boundaries if psp in passed in
