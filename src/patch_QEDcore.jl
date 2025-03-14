# TODO: write tests for this!

QEDcore.coordinate_name(::Type{Energy{IDX}}) where {IDX} = "energy_$IDX"

@generated function coordinate_symbol(
    coord::QEDcore.AbstractSingleParticleCoordinate{IDX},
) where {IDX}
    base_name = lowercase(string(nameof(coord)))
    coord_sym_local = base_name * "_$IDX"
    return quote
        Base.@_inline_meta
        (Symbol($coord_sym_local))
    end
end

@inline function coordinate_symbols(coord::QEDcore.AbstractUnivariateCoordinate)
    (coordinate_symbol(coord),)
end

coordinate_symbols(in_psl::TwoBodyRestSystem) = coordinate_symbols(in_psl.coord)
@inline coord_index(c_sym, coord_syms)::Int64 = findfirst(x -> x == c_sym, coord_syms)
