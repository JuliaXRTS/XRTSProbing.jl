function _build_coords(vp::XRTSProbing.VegasProposal, y::AbstractArray)
    coords = XRTSProbing._vegas_map(vp.vgrid, y)
    JAC = XRTSProbing._jac_vegas_map(vp.vgrid, y)

    return coords, JAC
end

function _generate_weights(
        rng::AbstractRNG,
        stp::AbstractProcessSetup,
        vp::XRTSProbing.VegasProposal,
        N::Int,
    )
    return _build_weights(stp, vp, rand(rng, N, ndims(vp)))
end

function _build_weights(
        stp::AbstractProcessSetup,
        vp::XRTSProbing.VegasProposal,
        rand_u::AbstractArray,
    )
    coords, jac = _build_coords(vp, rand_u)
    psps =
        PhaseSpacePoint.(
        process(vp),
        model(vp),
        phase_space_definition(vp),
        coords,
    )
    weights = @. compute(stp, psps) * jac

    return weights
end
