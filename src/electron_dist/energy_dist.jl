# TODO: make this a Single particle distribution (similar to the photon dists)

# source: https://en.wikipedia.org/wiki/Maxwell–Boltzmann_distribution#Distribution_for_the_momentum_vector
struct MaxellElectronEnergyDistribution{T, D} <: AbstractElectronEnergyDistribution
    dist::D
    temp::T

    function MaxellElectronEnergyDistribution(temp::T1) where {T <: Real, T1 <: Union{T, Quantity{T}}}
        temp_internal = _internalize_temperature(temp)
        a = sqrt(temp_internal) # mass = 1.0 for Electron
        dist = QEDevents.MaxwellBoltzmann(a)
        return new{T, typeof(dist)}(dist, temp_internal)
    end
end

temperature(en_dist::MaxellElectronEnergyDistribution) = en_dist.temp

function _energy_weight(en_dist::MaxellElectronEnergyDistribution, energy::T) where {T <: Real}
    rho = sqrt(energy^2 - one(energy))
    return pdf(en_dist.dist, rho)
end

energy_mean(en_dist::MaxellElectronEnergyDistribution{T}) where {T <: Real} = sqrt(mean(en_dist.dist)^2 + one(T))

# FIXME: QEDevents/src/patch_Distributions.jl:61, the theta needs to be an a
energy_width(en_dist::MaxellElectronEnergyDistribution{T}) where {T <: Real} = sqrt(var(en_dist.dist)^2 + one(T))
