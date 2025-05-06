
"""

    AbstractAxis

Abstract base type to describe the axis of a vector, e.g. the k-vector of a photon. Mostly used for multiple dispatch.
"""
abstract type AbstractAxis end
abstract type AbstractDefiniteAxis <: AbstractAxis end
abstract type AbstractIndefiniteAxis <: AbstractAxis end

struct XAxis <: AbstractDefiniteAxis end
struct YAxis <: AbstractDefiniteAxis end
struct ZAxis <: AbstractDefiniteAxis end

"""

    UniformAxis()

Helper type to indicate a random direction for the respective photon.
"""
struct UniformAxis <: AbstractIndefiniteAxis end

"""

    AbstractPhotonEnergyDistribution

Abstract base type for photon energy distributions. Subtypes of this type need to implement at least

* `energy_mean(dist::AbstractPhotonEnergyDistribution)`: return the mean energy modelled with `dist`
* `energy_width(dist::AbstractPhotonEnergyDistribution)`: return the energy width modelled with `dist`
* `k_vec_axis(dist::AbstractPhotonEnergyDistribution)`: return the direction (subtype of `AbstractAxis`) for the k-vector of the photons modelled with `dist`
* `_energyweight(dist::AbstractPhotonEnergyDistribution,omega)`: return the weight of `dist` for given photon energy `omega`
* `_randenergy(rng::AbstractRNG,dist::AbstractPhotonEnergyDistribution)`: return a random sample drawn from `dist`

Optionally, one may implement the following function to overwrite the generic implementation:

* `QEDevents._weight(dist::AbstractPhotonEnergyDistribution,mom)`: return the total weight for a given `mom` such that the sum over all possible moms results in unity.
* `QEDevents._randmom(rng::AbstractRNG,dist::AbstractPhotonEnergyDistribution)`: return a randomly distributed four-momentum according to `dist`

"""
abstract type AbstractPhotonEnergyDistribution{A<:AbstractAxis} <:
              SingleParticleDistribution end
@inline QEDevents._particle(d::AbstractPhotonEnergyDistribution) = Photon()
@inline QEDevents._particle_direction(d::AbstractPhotonEnergyDistribution) = Incoming()

# acutally, this needs to be different for different axes:
# - for axis::AbstractDefiniteAxis the k-vector axis weight is one
# - for axis::UniformAxis the k-vector axis weight is defined in such a way, that the
# integral over all possible momenta is one.
#
# However, since this adds a constant factor in the worst case, we will ignore this for
# now.
@inline QEDevents._weight(d::AbstractPhotonEnergyDistribution, omega) =
    _energyweight(d, omega)


function QEDevents._randmom(rng::AbstractRNG, d::AbstractPhotonEnergyDistribution{XAxis})
    energy = _randenergy(rng, d)
    return SFourMomentum(energy, energy, zero(energy), zero(energy))
end
function QEDevents._randmom(rng::AbstractRNG, d::AbstractPhotonEnergyDistribution{YAxis})
    energy = _randenergy(rng, d)
    return SFourMomentum(energy, zero(energy), energy, zero(energy))
end
function QEDevents._randmom(rng::AbstractRNG, d::AbstractPhotonEnergyDistribution{ZAxis})
    energy = _randenergy(rng, d)
    return SFourMomentum(energy, zero(energy), zero(energy), energy)
end
function QEDevents._randmom(
    rng::AbstractRNG,
    d::AbstractPhotonEnergyDistribution{UniformAxis},
)
    energy = _randenergy(rng, d)
    cth = 2 * rand(rng) - one(energy)
    sth = sqrt(one(energy) - cth^2)
    phi = 2 * pi * rand(rng)
    sph, cph = sin(phi), cos(phi)
    return SFourMomentum(energy, energy * sth * cph, energy * sth * sph, energy * cth)
end

# Gaussian distributed Photon energies

struct GaussianPhotonDist{T<:Real,D,A} <: AbstractPhotonEnergyDistribution{A}
    energy_mean::T
    energy_width::T
    dist::D
    k_vec_axis::A

    function GaussianPhotonDist(
        energy_mean::T,
        energy_width::T,
        k_vec_axis::A,
    ) where {T<:Real,A<:AbstractAxis}
        dist = truncated(Normal(energy_mean, energy_width), zero(T), _typed_inf(T))

        return new{T,typeof(dist),A}(energy_mean, energy_width, dist, k_vec_axis)
    end
end
@inline GaussianPhotonDist(energy_mean::T, energy_width::T) where {T<:Real} =
    GaussianPhotonDist(energy_mean, energy_width, UniformAxis())
energy_mean(d::GaussianPhotonDist) = d.energy_mean
energy_width(d::GaussianPhotonDist) = d.energy_width
k_vec_axis(d::GaussianPhotonDist) = d.k_vec_axis

@inline _energyweight(d::GaussianPhotonDist, omega::Real) = pdf(d.dist, omega)
@inline _randenergy(rng::AbstractRNG, d::GaussianPhotonDist) = rand(d.dist)


# uniform photon spectrum
#=
# add this after probing setup is merged
struct UniformPhotonDist{T<:Real,A,D} <: AbstractPhotonEnergyDistribution{A}
    max_energy::T
    k_vec_axis::A
    dist::D

    function UniformPhotonDist(
        max_energy::T,
        k_vec_axis::A,
    ) where {T<:Real,A<:AbstractAxis}

        dist = Uniform(zero(T),max_energy)
        return new{T,A,typeof(dist)}(max_energy,k_vec_axis,dist)

    end
end

k_vec_axis(d::UniformPhotonDist) = d.k_vec_axis
max_energy(d::UniformPhotonDist) = d.max_energy
@inline function _energyweight(d::UniformPhotonDist{T}, omega::S) where {T<:Real,S<:Real}
    zero(T)<= omega <= max_energy(omega) ? one(promote_type(T,S)) : zero(promote_type(T,S))
end
@inline _randenergy(rng::AbstractRNG, d::UniformPhotonDist) = rand(rng,d.dist)
=#

# TODO: implement other distributions
# - Poisson
# - Voigt
# - Lorentzian
# - data-driven (?)
