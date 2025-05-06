_typed_inf(::Type{T}) where {T <: Real} = typemax(T)
_typed_inf(::Type{Float64}) = Inf64
_typed_inf(::Type{Float32}) = Inf32
_typed_inf(::Type{Float16}) = Inf16
