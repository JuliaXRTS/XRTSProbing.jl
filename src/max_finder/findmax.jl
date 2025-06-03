function Base.findmax(
        rng::AbstractRNG,
        stp::AbstractProcessSetup,
        method::AbstractSampleBasedMaxFinder,
        sampler::ScatteringProcessDistribution,
    )

    # - build samples (momenta)
    # - build psps
    # - calculate weights==dcs
    # - perform maxfinder on weight array

    N = method.nsamples

    # TODO: is this general enough?
    samples, jac = _generate_coords(rng, sampler, N)

    # think about a three argument compute: compute(stp,sampler,samples)??
    weights = @. _compute(stp, samples) * jac

    return _findmax(method, weights)

end
