# implement rejection sampler here

struct EventGenerator{
    DCS<:DifferentialCrossSectionCached,
    PROP<:AbstractProposalDistribution,
    T<:Real,
} <: ScatteringProcessDistribution
    dcs::DCS
    proposal::PROP
    max_weight::T
end

QEDbase.process(eg::EventGenerator) = process(eg.dcs)
QEDbase.model(eg::EventGenerator) = model(eg.dcs)
QEDbase.phase_space_definition(eg::EventGenerator) = phase_space_definition(eg.dcs)

function generate_event(rng::AbstractRNG, eg::EventGenerator)
    rejected = true

    psp = PhaseSpacePoint(
        process(eg.dcs),
        model(eg.dcs),
        phase_space_definition(eg.dcs),
        Tuple(rand(SFourMomentum, number_incoming_particles(process(eg.dcs)))),
        Tuple(rand(SFourMomentum, number_outgoing_particles(process(eg.dcs)))),
    )
    w = 0.0

    while rejected
        rand_u = rand(rng, ndims(eg.proposal))
        rand_probability = rand(rng)

        coords, jac = _build_coords(eg.proposal, rand_u)

        psp = PhaseSpacePoint(
            process(eg.dcs),
            model(eg.dcs),
            phase_space_definition(eg.dcs),
            eg.dcs.in_coords,
            coords,
        )

        w = eg.dcs(psp) * jac

        if _rejection_filter(w, rand_probability, eg.max_weight)
            rejected = false
        end
    end

    return Event(psp, _residual_weight(w, eg.max_weight))
end

function _generate_weights_batch(rng::AbstractRNG, eg::EventGenerator, batch_size::Int)

    rand_us = rand(rng, batch_size, ndims(eg.proposal))
    rand_probabilies = rand(rng, batch_size)

    # TODO: consider using `sample(rng,dcs,sampler)` returning the psp and the weight
    coords, jac = _build_coords(eg.proposal, rand_us)
    psps =
        PhaseSpacePoint.(
            process(eg.proposal),
            model(eg.proposal),
            phase_space_definition(eg.proposal),
            Ref(eg.dcs.in_coords),
            coords,
        )
    weights = @. eg.dcs(psps) * jac

    mask = _rejection_filter.(weights, rand_probabilies, eg.max_weight)
    accepted_psps = _select_accepted(mask, psps)
    accepted_weights = _select_accepted(mask, weights)
    _update_residual_weight!(accepted_weights, eg.max_weight)

    return accepted_psps, accepted_weights
end

function generate_events(
    rng::AbstractRNG,
    eg::EventGenerator,
    n::Int,
    batch_size::Int = 1000,
)

    # TODO:
    # - consider redirect n=1 to the single event generation at compile time by using
    # Val(n)

    # TODO:
    # - put this into a function: alloc_event_list
    # - consider really allocating memory here, and not using sizehint!
    # - consider using shared memory for this to write the result in parallel
    accepted_weights = Float64[]
    psp_type = QEDevents._assemble_psp_type(
        process(eg.dcs),
        model(eg.dcs),
        phase_space_definition(eg.dcs),
        SFourMomentum,
    )
    accepted_psps = psp_type[]
    sizehint!(accepted_weights, n + batch_size - 1)
    sizehint!(accepted_psps, n + batch_size - 1)

    # TODO:
    # - this will be run in parallel
    nrun = 0
    while nrun <= n
        aaccepted_psps_batch, accepted_weights_batch =
            _generate_weights_batch(rng, eg, batch_size)

        # TODO: using scanreduce here to write in shared memory
        append!(accepted_weights, accepted_weights_batch)
        append!(accepted_psps, aaccepted_psps_batch)
        nrun += length(accepted_weights_batch)
    end

    # TODO:
    # - consider measuring the acceptence rate on the fly
    # - find out a way to opt this out (maybe similar to input validation)

    # TODO: is this final copy really necessary? (maybe doing everything on event level
    # is better?)
    return Event.(accepted_psps[1:n], accepted_weights[1:n])
end
