abstract type AbstractKinematicCuts{N,T} end
Base.broadcastable(c::AbstractKinematicCuts) = Ref(c)
degree_of_freedom(c::AbstractKinematicCuts{N}) where {N} = N

function _assert_correct_boundaries(::Tuple{}, ::Tuple{}) end

function _assert_correct_boundaries(
    low::Tuple{Vararg{T,N}},
    up::Tuple{Vararg{T,N}},
) where {T<:Real,N}
    first(low) <= first(up) || throw(
        ArgumentError(
            "lower boundary need to be smaller or equal to the respective upper boundary",
        ),
    )
    return _assert_correct_boundaries(low[2:end], up[2:end])
end

struct KinematicCuts{N,T} <: AbstractKinematicCuts{N,T}
    lower::NTuple{N,T}
    upper::NTuple{N,T}

    function KinematicCuts(low::NTuple{N,T}, up::NTuple{N,T}) where {N,T<:Real}
        _assert_correct_boundaries(low, up)

        return new{N,T}(low, up)
    end
end


Base.minimum(c::KinematicCuts) = c.lower
Base.maximum(c::KinematicCuts) = c.upper

Base.extrema(c::KinematicCuts) = (c.lower, c.upper)

function is_within_cuts(c::KinematicCuts{N}, d::NTuple{N}) where {N}
    minimum(c) .<= d .<= maximum(c)
end

function all_within_cuts(c::KinematicCuts{N}, d::NTuple{N}) where {N}
    all(is_within_cuts(c, d))
end
