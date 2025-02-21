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
