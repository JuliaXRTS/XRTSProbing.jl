
_zero_temperature_betabar(::Type{T}) where {T<:Real} = Inf
_zero_temperature_betabar(::Type{Float64}) = Inf64
_zero_temperature_betabar(::Type{Float32}) = Inf32
_zero_temperature_betabar(::Type{Float16}) = Inf16

_upper_bound_betabar(::Type{T}) where {T<:Real} = inv(eps(T))
