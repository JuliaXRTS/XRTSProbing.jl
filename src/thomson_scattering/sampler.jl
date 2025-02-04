
# TODO: introduce `method`, where inversion sampling could be another possibility
# (the current one could be called flat rejection)
struct SimpleThomsonSampler{T} <: ScatteringProcessDistribution
    # TODO: This will go to the PSL!
    omega::T
end

QEDbase.process(d::SimpleThomsonSampler) = ThomsonScattering()
QEDbase.model(d::SimpleThomsonSampler) = PerturbativeQED()
QEDbase.phase_space_definition(d::SimpleThomsonSampler) = ElabPhotonSphSystem()

function max_weight()
    return ALPHA_SQUARE / ELECTRONMASS^2
end

# trivial rejection sampler -> consider implementing inversion sampler!
function _TS_generate_coord_event_serial(rng::AbstractRNG)
    maximum_weight = max_weight()
    while true
        cth_trail = 2 * rand(rng) - 1

        weight = _TS_diffCS_pol_spin_sum(cth_trail) / maximum_weight

        if weight >= rand(rng)
            phi_trail = 2 * pi * rand(rng)
            return (cth = cth_trail, phi = phi_trail)
        end
    end
end

function QEDevents._randmom(rng::AbstractRNG, d::SimpleThomsonSampler)
    coord_event = _TS_generate_coord_event_serial(rng)
    om = d.omega

    in_moms = QEDbase._generate_incoming_momenta(
        process(d),
        model(d),
        phase_space_definition(d),
        (om,),
    )

    out_moms = QEDbase._generate_outgoing_momenta(
        process(d),
        model(d),
        phase_space_definition(d),
        (om,),
        (coord_event.cth, coord_event.phi),
    )
    return (in_moms, out_moms)
end

function QEDevents._weight(d::SimpleThomsonSampler, psp::PhaseSpacePoint)
    return unsafe_differential_cross_section(psp)
end
