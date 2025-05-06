function insert_by_name(syms, rem::Tuple{Vararg{T}}; kwargs...) where {T <: Real}
    N = length(syms)
    k = keys(kwargs)
    l = length(k)
    k_pos = ntuple(i -> coord_index(k[i], syms), Val(l))
    reduced_syms = TupleTools.deleteat(syms, k_pos)

    return ntuple(Val(N)) do i
        s = syms[i]
        if s in k
            return kwargs[s]
        else
            idx = coord_index(s, reduced_syms)
            return rem[idx]
        end
    end
end
