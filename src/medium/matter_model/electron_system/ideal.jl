# TODO:
# - move Lindhard
# - remove folder structure_factor

struct IdealElectronSystem{TM, T} <: AbstractProperElectronSystem
    electron_density::T
    temperature::T

    function IdealElectronSystem{TM}(
            electron_density::T1,
            temperature::T2,
        ) where {
            T <: Real,
            T1 <: Union{T, Quantity{T}},
            T2 <: Union{T, Quantity{T}},
            TM <: AbstractTemperatureMode,
        }

        ne_internal = _internalize_density(electron_density)
        temp_internal = _internalize_temperature(temperature)

        return new{TM, T}(ne_internal, temp_internal)
    end
end

function IdealElectronSystem{ZeroTemperature}(
        electron_density::T1,
    ) where {T <: Real, T1 <: Union{T, Quantity{T}}}
    return IdealElectronSystem{ZeroTemperature}(electron_density, zero(T))
end

function IdealElectronSystem(
        electron_density::T1,
        temperature::T2,
    ) where {T <: Real, T1 <: Union{T, Quantity{T}}, T2 <: Union{T, Quantity{T}}}
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature)
end

@inline temperature(elsys::IdealElectronSystem{ZeroTemperature, T}) where {T <: Real} = zero(T)
@inline temperature(elsys::IdealElectronSystem{FiniteTemperature, T}) where {T <: Real} =
    elsys.temperature
@inline electron_density(elsys::IdealElectronSystem) = elsys.electron_density


function imag_dynamic_response(
        elsys::IdealElectronSystem{ZeroTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    return imag_lindhard_zero_temperature(ombar, qbar)
end

function real_dynamic_response(
        elsys::IdealElectronSystem{ZeroTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    return real_lindhard_zero_temperature(ombar, qbar)
end

function imag_dynamic_response(
        elsys::IdealElectronSystem{FiniteTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)

    if ombar <= zero(ombar)
        if qbar < zero(qbar)
            return imag_lindhard_nonzero_temperature(-ombar, -qbar, betabar(elsys))
        else
            return -imag_lindhard_nonzero_temperature(-ombar, qbar, betabar(elsys))
        end
    end

    return imag_lindhard_nonzero_temperature(ombar, qbar, betabar(elsys))
end

function real_dynamic_response(
        elsys::IdealElectronSystem{FiniteTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    return real_lindhard_nonzero_temperature(ombar, qbar, betabar(elsys))
end

# zero temperature

function _imag_lindhardDSF_zeroT_minus(omb::T, qb::T) where {T <: Real}
    num = _nu_minus(omb, qb)

    if abs(num) >= one(T)
        return zero(T)
    end
    return one(omb) - num^2

end

function _imag_lindhardDSF_zeroT_plus(omb::T, qb::T) where {T <: Real}
    nup = _nu_plus(omb, qb)
    if abs(nup) >= one(T)
        return zero(T)
    end
    return one(omb) - nup^2

end

function imag_lindhard_zero_temperature(ombar::T, qbar::T) where {T <: Real}
    return -pi / (4 * qbar) *
        (_imag_lindhardDSF_zeroT_minus(ombar, qbar) - _imag_lindhardDSF_zeroT_plus(ombar, qbar))
end

function real_lindhard_zero_temperature(ombar::T, qbar::T) where {T <: Real}
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    term1 = (1 - num^2) / 4 * _stable_log_term(num)
    term2 = (1 - nup^2) / 4 * _stable_log_term(nup)

    return -(qbar / 2 - term1 + term2) / qbar
end

# finite temperature

function _integrand_real_Lindhard_finT(omb, qb, x, bbar)
    num = QEDprobing._nu_minus(omb, qb)
    nup = QEDprobing._nu_plus(omb, qb)
    term1 = log(abs((x - num) / (x + num)))
    term2 = log(abs((x - nup) / (x + nup)))
    return _F(x, bbar) * (term1 - term2) / 2
end

function real_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T) where {T <: Real}

    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)

    res, _ = quadgk(
        x -> _integrand_real_Lindhard_finT(ombar, qbar, x, bbar),
        zero(num),
        max(num, zero(num)),
        nup,
        Inf,
    )

    return -res / (qbar)
end

function _imag_lindhard_large_om(omb, qb, bbar)
    num = QEDprobing._nu_minus(omb, qb)
    prefac = -pi / (4 * qb * bbar)
    mubar = QEDprobing._chemical_potential_normalized(bbar)
    return prefac * exp(-bbar * (num^2 - mubar))
end

function _om_crit_func(qbar, bbar)
    e = 1.0e-8
    mubar = QEDprobing._chemical_potential_normalized(bbar)
    log_arg = 4 * e / pi * qbar * bbar
    r = abs(mubar - inv(bbar) * log(log_arg))
    return 2 * qbar * sqrt(r) + qbar^2
end

function _imag_lindhard_finT_unstable(ombar, qbar, bbar)
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = log1pexp(bbar * (num^2 - mubar))
    term2 = log1pexp(bbar * (nup^2 - mubar))

    return -pi / (4 * qbar) * (ombar + (term1 - term2) / bbar)
end


function imag_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T) where {T <: Real}
    if bbar < 5.0 * one(T) &&
            zero(qbar) < abs(qbar) <= 2.0e-3 &&
            abs(ombar) >= _om_crit_func(qbar, bbar)
        return _imag_lindhard_large_om(ombar, qbar, bbar)
    else
        return _imag_lindhard_finT_unstable(ombar, qbar, bbar)
    end
end
