struct NaiveMaxFinder <: AbstractSampleBasedMaxFinder
    nsamples::Int
end

function _findmax(method::NaiveMaxFinder, weights::AbstractVector)
    return maximum(weights)
end
