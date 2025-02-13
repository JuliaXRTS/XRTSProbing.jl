

struct QuantileReductionMethod <: AbstractSampleBasedMaxFinder
    p::Float64
    nsamples::Int
end

function _findmax(method::QuantileReductionMethod, weights::AbstractVector)
    sorted_weights = sort(weights)
    sum_weight_quantile = method.p * sum(sorted_weights)
    s = 0.0
    out_weight = 0.0
    for (i, weight) in enumerate(sorted_weights)
        s += weight
        if s > sum_weight_quantile
            out_weight = sorted_weights[i-1]
        end
    end

    return out_weight
end

function _findmax1(method::QuantileReductionMethod, weights::AbstractVector)
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
end
