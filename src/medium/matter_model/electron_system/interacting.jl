struct InteractingElectronSystem{PES<:AbstractProperElectronSystem,S<:AbstractScreening} <:
       AbstractInteractingElectronSystem
    proper_elec_system::PES
    screening::S

    function InteractingElectronSystem(
        pelsys::PES,
        scr::S,
    ) where {PES<:AbstractProperElectronSystem,S<:AbstractScreening}

        return new{PES,S}(pelsys, scr)
    end

    function InteractingElectronSystem{PES}(
        electron_density::T1,
        temp::T2,
        scr::S,
    ) where {
        T<:Real,
        T1<:Union{T,Quantity},
        T2<:Union{T,Quantity},
        PES<:AbstractProperElectronSystem,
        S<:AbstractScreening,
    }
        proper_elec_system = PES(electron_density, temp)
        return new{PES,S}(proper_elec_system, scr)
    end
end

# accessors
@inline proper_electron_system(elsys::InteractingElectronSystem) = elsys.proper_elec_system
@inline screening(elsys::InteractingElectronSystem) = elsys.screening
@inline temperature(elsys::InteractingElectronSystem) =
    temperature(proper_electron_system(elsys))
@inline electron_density(elsys::InteractingElectronSystem) =
    electron_density(proper_electron_system(elsys))

# delegations
@inline function local_effective_potential(
    iesys::AbstractInteractingElectronSystem,
    om_q::NTuple{2,T},
) where {T<:Real}
    return local_effective_potential(screening(iesys), om_q)
end

@inline function pseudo_potential(
    iesys::AbstractInteractingElectronSystem,
    om_q::NTuple{2,T},
) where {T<:Real}
    return pseudo_potential(screening(iesys), om_q)
end

@inline function dielectric_function(
    iesys::AbstractInteractingElectronSystem,
    om_q::NTuple{2,T},
) where {T<:Real}
    return dielectric_function(proper_electron_system(iesys), screening(iesys), om_q)
end

# response functions

function dynamic_response(
    iesys::AbstractInteractingElectronSystem,
    om_q::NTuple{2,T},
) where {T<:Real}
    rf = dynamic_response(proper_electron_system(iesys), om_q)
    lep = local_effective_potential(screening(iesys), om_q)
    eps = _dielectric_function(rf, lep)
    return rf / eps
end

function real_dynamic_response(
    iesys::AbstractInteractingElectronSystem,
    om_q::NTuple{2,T},
) where {T<:Real}
    rf = dynamic_response(iesys, om_q)
    return real(rf)
end

function imag_dynamic_response(
    iesys::AbstractInteractingElectronSystem,
    om_q::NTuple{2,T},
) where {T<:Real}
    rf = dynamic_response(iesys, om_q)
    return imag(rf)
end
