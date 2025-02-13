
function findmax(
    rng::AbstractRNG,
    dcs::DifferentialCrossSectionCached,
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

    weights = @. dcs(samples) * jac

    return _findmax(method, weights)

end
