
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

    AllAxis()

Helper type to indicate a random direction for the respective photon.
"""
struct AllAxis <: AbstractIndefiniteAxis end

abstract type AbstractPhotonDistribution <: SingleParticleDistribution end
@inline QEDevents._particle(d::AbstractPhotonDistribution) = Photon()
@inline QEDevents._particle_direction(d::AbstractPhotonDistribution) = Incoming()

# Gaussian distributed Photon energies

struct GaussianPhotonDist{T<:Real,D,A} <: AbstractPhotonDistribution
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
    GaussianPhotonDist(energy_mean, energy_width, AllAxis())
energy_mean(d::GaussianPhotonDist) = d.energy_mean
energy_width(d::GaussianPhotonDist) = d.energy_width
k_vec_axis(d::GaussianPhotonDist) = d.k_vec_axis

@inline QEDevents._weight(d::GaussianPhotonDist, x) = pdf(d.dist, x)

function QEDevents._randmom(rng::AbstractRNG, d::GaussianPhotonDist{T,D,XAxis}) where {T,D}
    energy = rand(rng, d.dist)
    return SFourMomentum(energy, energy, zero(T), zero(T))
end
function QEDevents._randmom(rng::AbstractRNG, d::GaussianPhotonDist{T,D,YAxis}) where {T,D}
    energy = rand(rng, d.dist)
    return SFourMomentum(energy, zero(T), energy, zero(T))
end
function QEDevents._randmom(rng::AbstractRNG, d::GaussianPhotonDist{T,D,ZAxis}) where {T,D}
    energy = rand(rng, d.dist)
    return SFourMomentum(energy, zero(T), zero(T), energy)
end
function QEDevents._randmom(
    rng::AbstractRNG,
    d::GaussianPhotonDist{T,D,AllAxis},
) where {T,D}
    energy = rand(rng, d.dist)
    cth = 2 * rand(rng) - one(T)
    #sth = sqrt(one(T)-cth^2)
    sth = sin(acos(cth))
    phi = 2 * pi * rand(rng)
    sph, cph = sin(phi), cos(phi)
    return SFourMomentum(energy, energy * sth * cph, energy * sth * sph, energy * cth)
end

# TODO: implement other distributions
# - Poisson
# - Voigt
# - Lorentzian
# - data-driven (?)
