# implement rejection sampler here

struct EventGenerator{
        STP <: AbstractProcessSetup,
        PROP <: AbstractProposalDistribution,
        T <: Real,
    } <: ScatteringProcessDistribution
    stp::STP
    proposal::PROP
    max_weight::T
end

setup(eg::EventGenerator) = eg.stp
proposal(eg::EventGenerator) = eg.proposal
max_weight(eg::EventGenerator) = eg.max_weight
QEDbase.process(eg::EventGenerator) = process(setup(eg))
QEDbase.model(eg::EventGenerator) = model(setup(eg))
QEDbase.phase_space_layout(eg::EventGenerator) = phase_space_layout(setup(eg))

function generate_event(rng::AbstractRNG, eg::EventGenerator)
    rejected = true

    # build a psp just to allocate the memory
    psp = PhaseSpacePoint(
        process(eg),
        model(eg),
        phase_space_layout(eg),
        Tuple(rand(SFourMomentum, number_incoming_particles(process(eg)))),
        Tuple(rand(SFourMomentum, number_outgoing_particles(process(eg)))),
    )
    w = 0.0

    while rejected
        rand_u = rand(rng, ndims(proposal(eg)))
        rand_probability = rand(rng)

        coords, jac = _build_coords(proposal(eg), rand_u)

        psp = _build_psp(proposal(eg), coords)
        w = compute(setup(eg), psp) * jac

        if _rejection_filter(w, rand_probability, max_weight(eg))
            rejected = false
        end
    end

    return Event(psp, _residual_weight(w, max_weight(eg)))
end

function _generate_weights_batch(rng::AbstractRNG, eg::EventGenerator, batch_size::Int)


    # generation of rnd args
    ## generate vegas coords
    rand_us = rand(rng, batch_size, ndims(proposal(eg)))
    rand_probabilies = rand(rng, batch_size)

    # TODO: consider using `sample(rng,dcs,sampler)` returning the psp and the weight
    #transform vegas coords into setup coords
    coords, jac = _build_coords(proposal(eg), rand_us)

    # build psp from setup coords
    psps = _build_psp.(proposal(eg), coords)

    # rng vals
    # build full weight from setup coords and proposal jac
    weights = @. compute(setup(eg), psps) * jac

    # filter mask
    mask = _rejection_filter.(weights, rand_probabilies, max_weight(eg))

    # arg selection
    accepted_psps = _select_accepted(mask, psps)

    # value selection
    accepted_weights = _select_accepted(mask, weights)

    # value update
    _update_residual_weight!(accepted_weights, max_weight(eg))

    # return
    return accepted_psps, accepted_weights
end

function generate_events(
        rng::AbstractRNG,
        eg::EventGenerator,
        n::Int,
        batch_size::Int = min(n, 1000),
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
        process(eg),
        model(eg),
        phase_space_layout(eg),
        SFourMomentum{Float64},
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
