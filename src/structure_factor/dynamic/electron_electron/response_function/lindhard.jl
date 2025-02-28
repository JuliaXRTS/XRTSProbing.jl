
abstract type AbstractLindhardDRF{TEMPMODE<:AbstractTemperatureMode,T<:Real} <:
              AbstractDynamicResponseFunction{TEMPMODE,T} end


"""

    LindhardDRF(fermi_wave_vector_mag::Real)
    LindhardDRF(temp::Real,fermi_wave_vector_mag::Real)

Return the Lindhard dynamic response function.

!!! note "Normalization"

    The current implementation drops the constant normalization, i.e. a call returns \$\\chi_0(\\omega,q)/N(0)\$.

"""
struct LindhardDRF{TMODE,T<:Real} <: AbstractLindhardDRF{TMODE,T}
    betabar::T

    function LindhardDRF(::Type{T}) where {T<:Real}
        return new{ZeroTemperature,T}(_zero_temperature_betabar(T))
    end

    function LindhardDRF(betabar::T) where {T<:Real}
        return new{FiniteTemperature,T}(min(betabar, _upper_bound_betabar(T)))
    end
end

LindhardDRF() = LindhardDRF(Float64)

betabar(rf::LindhardDRF) = getfield(rf, :betabar)

# zero temperature

function _imag_lindhardDSF_zeroT_minus(omb::T, qb::T) where {T<:Real}
    num = _nu_minus(omb, qb)

    if abs(num) >= one(T)
        return zero(T)
    end
    return one(omb) - num^2

end

function _imag_lindhardDSF_zeroT_plus(omb::T, qb::T) where {T<:Real}
    nup = _nu_plus(omb, qb)
    if abs(nup) >= one(T)
        return zero(T)
    end
    return one(omb) - nup^2

end

function _compute_imag(
    rf::LindhardDRF{ZeroTemperature},
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    ombar, qbar = ombarqbar
    -pi / (2 * qbar) *
    (_imag_lindhardDSF_zeroT_minus(ombar, qbar) - _imag_lindhardDSF_zeroT_plus(ombar, qbar))
end

function _compute_real(
    rf::LindhardDRF{ZeroTemperature},
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    ombar, qbar = ombarqbar
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    term1 = (1 - num^2) / (2 * qbar) * log(abs((num + 1) / (num - 1)))
    term2 = (1 - nup^2) / (2 * qbar) * log(abs((nup + 1) / (nup - 1)))

    return -(1 - term1 + term2)
end


# finite temperature

function _integrand_real_Lindhard_finT(omb, qb, x, bbar)
    num = QEDprobing._nu_minus(omb, qb)
    nup = QEDprobing._nu_plus(omb, qb)
    term1 = log(abs((x - num) / (x + num)))
    term2 = log(abs((x - nup) / (x + nup)))
    return _F(x, bbar) * (term1 - term2)
end

function _compute_real(
    rf::LindhardDRF{FiniteTemperature},
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    ombar, qbar = ombarqbar
    bbar = betabar(rf)

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

function _compute_imag(
    rf::LindhardDRF{FiniteTemperature},
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    ombar, qbar = ombarqbar
    bbar = betabar(rf)
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = log1pexp(bbar * (num^2 - mubar))
    term2 = log1pexp(bbar * (nup^2 - mubar))

    return -pi / (2 * qbar) * (ombar + (term1 - term2) / bbar)
end
