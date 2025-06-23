struct Event{PSP <: AbstractPhaseSpacePoint, T <: Number}
    psp::PSP
    weight::T
end

QEDcore.process(e::Event) = process(e.psp)
QEDcore.model(e::Event) = model(e.psp)
QEDcore.phase_space_layout(e::Event) = phase_space_definition(e.psp)

phase_space_point(e::Event) = e.psp
weight(e::Event) = e.weight

### I/O

##### SAVING #####

function save(
        event_list::AbstractVector{E},
        path::String = "data.h5"
    ) where {
        E <: Event,
    }
    moms = @. momenta(phase_space_point(event_list))
    ws = weight.(event_list)

    first_event = first(event_list)
    h5open(path, "w") do file
        file["momenta"] = moms
        file["weights"] = ws
    end
    return nothing
end

function _get_direction(mom::AbstractFourMomentum)
    return (getX(mom), getY(mom), getZ(mom)) ./ getMag(mom)
end

function save_photon_en_mag(
        event_list::AbstractVector{E},
        path::String = "photon_en_mag.h5";
        keV = false,
    ) where {
        E <: Event,
    }
    photon_moms = @. momentum(phase_space_point(event_list), Outgoing(), Photon())
    om_prime = getE.(photon_moms)
    mom_dir = _get_direction.(photon_moms)

    scale = keV ? (ELECTRONMASS / 1.0e3) : one(Float64)
    h5open(path, "w") do file
        file["omega_prime"] = scale * om_prime
        file["direction"] = mom_dir
    end
    return nothing

end

# helper to convert
function _to_moms_tuple(nt::NamedTuple)
    moms_nt = values(nt)
    return moms = @. SFourMomentum(values(moms_nt))
end

function load(
        filepath::String,
        proc,
        model,
        psl
    )
    moms, ws = h5open(filepath, "r") do file
        moms = read(file, "momenta")
        ws = read(file, "weights")

        (_to_moms_tuple.(moms), ws)
    end


    psps = PhaseSpacePoint.(proc, model, psl, moms)
    event_list = Event.(psps, ws)
    return event_list
end
