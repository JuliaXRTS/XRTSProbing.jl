struct QuantileReductionMethod <: AbstractSampleBasedMaxFinder
    p::Float64
    nsamples::Int
end

function _findmax(method::QuantileReductionMethod, weights::AbstractVector)
    sorted_weights = sort(weights)
    sum_weight = sum(sorted_weights)
    sum_weight_quantile = method.p * sum_weight

    s = sum_weight
    for weight in sorted_weights
        s -= weight
        if s <= sum_weight_quantile
            return weight
        end
    end
    return sorted_weights[end]
end
