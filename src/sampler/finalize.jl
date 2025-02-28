function _residual_weight(accepted_weight::Real, maximum_weight::Real)
    return max(accepted_weight / maximum_weight, one(maximum_weight))
end

function _update_residual_weight!(
    accepted_weights::AbstractVector{T},
    maximum_weight::Real,
) where {T<:Real}
    for i in eachindex(accepted_weights)
        accepted_weights[i] = max(accepted_weights[i] / maximum_weight, one(maximum_weight))
    end
    return accepted_weights
end

# TODO: test if this is allocating
function _update_residual_weight!(
    accepted_events::AbstractVector{T},
    maximum_weight::Real,
) where {T<:Event}
    for i in eachindex(accepted_events)
        residual_weight = _residual_weight(accepted_events[i], maximum_weight)
        accepted_events[i] = Event(accepted_events[i].psp, residual_weight)
    end
    return accepted_events
end
