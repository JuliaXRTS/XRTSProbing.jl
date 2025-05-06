struct Vacuum <: AbstractMatter end

function dynamic_structure_factor(esys::Vacuum, om_q::NTuple{2, T}) where {T <: Real}
    return one(T)
end
