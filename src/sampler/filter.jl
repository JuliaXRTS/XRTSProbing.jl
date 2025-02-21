@inline function _rejection_filter(trial_weight::Real, prob::Real, w_max::Real)
    return trial_weight >= prob * w_max
end

@inline function _rejection_filter(
    trail_weights::AbstractVector,
    probs::AbstractVector,
    w_max::Real,
)
    # TODO: avoid such functions -> not idiomatic Julia
    return _rejection_filter.(trail_weights, probs, w_max)
end
