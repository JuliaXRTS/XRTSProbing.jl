abstract type AbstractResponseFunctionBasedEEDSF <: AbstractElectronElectronDSF end
"""

    response_function(dsf::AbstractResponseFunctionBasedEEDSF)

Interface function: return associated response function, i.e. subtype of `AbstractElectronElectronResponseFunction`.
"""
function response_function end

temperature(lh::AbstractResponseFunctionBasedEEDSF) = temperature(response_function(lh))

abstract type AbstractDynamicResponseFunction{TEMPMODE<:AbstractTemperatureMode,T<:Real} end

"""
    temperature(::AbstractDynamicResponseFunction)

Return the temperature value in associated with the given dynamic response function in units of the Fermi energy: \$k_BT/E_F\$

!!! note "zero temperature"

    For `ZeroTemperature` mode, this returns exactly zero, for `FiniteTemerature` mode, the temperature value
    is bounded from below by `eps()` to avoid `NaN` and `Inf`.

"""
@inline function temperature(
    ::AbstractDynamicResponseFunction{ZeroTemperature,T},
) where {T<:Real}
    return zero(T)
end

@inline function temperature(
    rf::AbstractDynamicResponseFunction{FiniteTemperature,T},
) where {T<:Real}
    return inv(betabar(rf))
end

"""

    betabar(::AbstractDynamicResponseFunction)

Interface function: return the value of \$\\bar\\beta:=E_F/k_BT\$, where \$T\$ denotes the temperature.

!!! note "zero temperature"

    For `ZeroTemperature` mode, this returns exactly `Inf`; for `FiniteTemerature` mode, the `betabar` value
    is bounded from below by `inv(eps())` to avoid `NaN` and `Inf`.

"""
function betabar(::AbstractDynamicResponseFunction) end

@inline function betabar(
    ::AbstractDynamicResponseFunction{ZeroTemperature,T},
) where {T<:Real}
    return _zero_temperature_betabar(T)
end

"""

    _compute_real(rf::AbstractDynamicResponseFunction,q_omega::Tuple)

Interface function: return real part of response function associated with `rf`.

"""
function _compute_real end

"""

    _compute_imag(rf::AbstractDynamicResponseFunction,q_omega::Tuple)

Interface function: return imaginary part of response function associated with `rf`.

"""
function _compute_imag end

@inline function _compute(
    rf::AbstractDynamicResponseFunction,
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    return _compute_real(rf, ombarqbar) + 1im * _compute_imag(rf, ombarqbar)
end

@inline function (rf::AbstractDynamicResponseFunction)(
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    return _compute(rf, ombarqbar)
end

@inline function Base.imag(
    rf::AbstractDynamicResponseFunction,
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    return _compute_imag(rf, ombarqbar)
end

@inline function Base.real(
    rf::AbstractDynamicResponseFunction,
    ombarqbar::NTuple{2,T},
) where {T<:Real}
    return _compute_real(rf, ombarqbar)
end

#abstract type AbstractPseudoPotentialBasedDRF{TEMPMODE<:AbstractTemperatureMode,T<:Real} <: AbstractDynamicResponseFunction{TEMPMODE,T} end

#=
# concrete implementation
struct ResponseFunctionBasedDSF{RF<:AbstractDynamicResponseFunction} <: AbstractResponseFunctionBasedEEDSF
    rf::RF
end

response_function(dsf::ResponseFunctionBasedDSF) = dsf.rf

# Lindhard dynamic structure factor

const LindhardDSF{T} = ResponseFunctionBasedDSF{RF} where {T<:Real,RF<:AbstractLindhardDRF{T}}
function LindhardDSF(
    temp::T,
    pol_deg::T,
    n_density::T
) where {T<:Real}
    rf = LindhardDRF(temp,pol_deg,n_density)
    return ResponseFunctionBasedDSF(rf)
end
=#
