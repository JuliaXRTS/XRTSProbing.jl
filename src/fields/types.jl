
abstract type AbstractSpectrumBasedField <: QEDfields.AbstractBackgroundField end
struct DistributionBasedField{D<:AbstractPhotonDistribution} <: AbstractSpectrumBasedField
    dist::D
end

energy_spectrum(f::DistributionBasedField, energy::T) where {T<:Real} =
    QEDevents._weight(f.dist, energy)
