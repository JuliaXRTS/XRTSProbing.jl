
function QEDevents._randmom(rng::AbstractRNG, vp::VegasProposal)
    D = ndims(vp)
    y = rand(rng, D)
    coords = _vegas_map(vp.vgrid, y)

    in_moms = QEDbase._generate_incoming_momenta(
        process(vp),
        model(vp),
        phase_space_definition(vp),
        vp.dcs.in_coords,
    )

    out_moms = QEDbase._generate_outgoing_momenta(
        process(vp),
        model(vp),
        phase_space_definition(vp),
        vp.dcs.in_coords,
        coords,
    )

    return in_moms, out_moms
end

function generate_events(rng::AbstractRNG, vp::VegasProposal)
    D = ndims(vp)
    y = rand(rng, D)
    return _generate_events(y, vp)
end

@inline function _generate_events(rand_u, vp::VegasProposal)
    coords = _vegas_map(vp.vgrid, rand_u)
    psp = PhaseSpacePoint(
        process(vp),
        model(vp),
        phase_space_definition(vp),
        vp.dcs.in_coords,
        coords,
    )

    JAC = _jac_vegas_map(vp.vgrid, rand_u)
    F = vp.dcs(psp)
    weight = F * JAC

    return Event(psp, weight)
end

function generate_events(rng::AbstractRNG, vp::VegasProposal, n::Int)
    D = ndims(vp)
    y = rand(rng, n, D)
    coords = _vegas_map(vp.vgrid, y)
    JAC = _jac_vegas_map(vp.vgrid, y)

    dest = Vector{Event{eltype(vp),Float64}}(undef, n)

    @inbounds for i in eachindex(dest)
        dest[i] = _build_event(vp.dcs, coords[i], JAC[i])
    end

    return dest
end
