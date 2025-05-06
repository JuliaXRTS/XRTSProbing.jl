# local effective potential


@inline function _local_effective_potential(vp, lfc)
    return vp * (one(lfc) - lfc)
end

function local_effective_potential(
        scr::AbstractScreening,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    vp = pseudo_potential(scr, om_q)
    lfc = local_field_correction(scr, om_q)

    return _local_effective_potential(vp, lfc)
end

"""

    dielectric_function(::AbstractIdealElectronSystem,::AbstractScreening,om_q)

Return the value of the dielectric function for the given proper response function and given local effective potential.
"""
function _dielectric_function(rf::Number, lep::Number)
    temp = lep * rf
    return one(temp) - temp
end

"""

    dielectric_function(::AbstractIdealElectronSystem,::AbstractScreening,om_q)

Return the value of the dielectric function for the given electron system and given screening.
"""
function dielectric_function(
        esys::AbstractProperElectronSystem,
        scr::AbstractScreening,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    rf = dynamic_response(esys, om_q)
    lep = local_effective_potential(scr, om_q)
    return _dielectric_function(rf, lep)
end

# local field correction

struct NoLocalFieldCorrection <: AbstractLocalFieldCorrection end

function _compute(::NoLocalFieldCorrection, om_q::NTuple{2, T}) where {T <: Real}
    return zero(T)
end

# pseudo potential

abstract type AbstractPseudoPotential end

struct CoulombPseudoPotential <: AbstractPseudoPotential end

function _compute(::CoulombPseudoPotential, om_q::NTuple{2, T}) where {T <: Real}
    om, q = om_q
    return ELEMENTARY_CHARGE_SQUARED / q^2
end

### screening

# no screening

struct NoScreening <: AbstractScreening end

dielectric_function(
    esys::AbstractProperElectronSystem,
    scr::NoScreening,
    om_q::NTuple{2, T},
) where {T <: Real} = one(T)

local_effective_potential(src::NoScreening, om_q::NTuple{2, T}) where {T <: Real} = zero(T)
@inline pseudo_potential(scr::NoScreening, om_q::NTuple{2, T}) where {T <: Real} = zero(T)
@inline local_field_correction(scr::NoScreening, om_q::NTuple{2, T}) where {T <: Real} =
    zero(T)

# screeing with pseudo potential and local field correction

struct Screening{PP <: AbstractPseudoPotential, LFC <: AbstractLocalFieldCorrection} <:
    AbstractScreening
    pseudo_potential::PP
    lfc::LFC
end

Screening() = Screening(CoulombPseudoPotential(), NoLocalFieldCorrection())

@inline pseudo_potential(scr::Screening, om_q) = _compute(scr.pseudo_potential, om_q)
@inline local_field_correction(scr::Screening, om_q) = _compute(scr.lfc, om_q)
