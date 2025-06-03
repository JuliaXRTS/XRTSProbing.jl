function _build_coords(vp::QEDprobing.VegasProposal, y::AbstractArray)
    coords = QEDprobing._vegas_map(vp.vgrid, y)
    JAC = QEDprobing._jac_vegas_map(vp.vgrid, y)

    return coords, JAC
end

function _generate_weights(
        rng::AbstractRNG,
        stp::AbstractProcessSetup,
        vp::QEDprobing.VegasProposal,
        N::Int,
    )
    return _build_weights(stp, vp, rand(rng, N, ndims(vp)))
end

function _build_weights(
        stp::AbstractProcessSetup,
        vp::QEDprobing.VegasProposal,
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
