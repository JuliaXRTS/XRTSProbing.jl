abstract type AbstractSpectrumBasedField <: QEDfields.AbstractBackgroundField end
struct DistributionBasedField{D <: AbstractPhotonEnergyDistribution} <:
    AbstractSpectrumBasedField
    dist::D
end

energy_spectrum(f::DistributionBasedField, energy::T) where {T <: Real} =
    QEDevents._weight(f.dist, energy)

energy_mean(f::DistributionBasedField) = energy_mean(f.dist)
energy_width(f::DistributionBasedField) = energy_width(f.dist)
