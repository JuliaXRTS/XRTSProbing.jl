abstract type AbstractElectronDistribution end
abstract type AbstractElectronEnergyDistribution end


struct EnergyBasedElectronDistribution{D <: AbstractElectronEnergyDistribution} <: AbstractElectronDistribution
    en_dist::D
end

temperature(dist::EnergyBasedElectronDistribution) = temperature(dist.en_dist)
energy_spectrum(dist::EnergyBasedElectronDistribution, energy::T) where {T <: Real} =
    _energy_weight(dist.en_dist, energy)

energy_mean(dist::EnergyBasedElectronDistribution) = energy_mean(dist.en_dist)
energy_width(dist::EnergyBasedElectronDistribution) = energy_width(dist.en_dist)
