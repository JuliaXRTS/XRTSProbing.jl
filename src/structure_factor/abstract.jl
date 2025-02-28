abstract type AbstractStructureFactor end

"""

    _compute(sf::AbstractStructureFactor,input)
    _compute(sf::AbstractResponseFunction,input)

Interface function: return structure factor associated with `sf`. With `T<:Real`, the argument
`input` can be `q::T` for static or `q_omega::NTuple{2,T}` for dynamic structure factors.

"""
function _compute end

abstract type AbstractDynamicStructureFactor <: AbstractStructureFactor end

function (dsf::AbstractDynamicStructureFactor)(q_omega::NTuple{2,T}) where {T<:Real}
    return _compute(dsf, q_omega)
end

abstract type AbstractSingleInteractionDSF <: AbstractDynamicStructureFactor end
abstract type AbstractElectronElectronDSF <: AbstractSingleInteractionDSF end
abstract type AbstractCompoundDSF <: AbstractDynamicStructureFactor end

abstract type AbstractStaticStructureFactor <: AbstractStructureFactor end
abstract type AbstractSingleInteractionSSF <: AbstractStaticStructureFactor end
abstract type AbstractCompoundSSF <: AbstractStaticStructureFactor end


### temperature dependence

# Temperature modes
abstract type AbstractTemperatureMode end
struct FiniteTemperature <: AbstractTemperatureMode end
struct ZeroTemperature <: AbstractTemperatureMode end
