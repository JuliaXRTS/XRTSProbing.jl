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
    return (coordinate_symbol(coord),)
end

coordinate_symbols(in_psl::TwoBodyRestSystem) = coordinate_symbols(in_psl.coord)
@inline coord_index(c_sym, coord_syms)::Int64 = findfirst(x -> x == c_sym, coord_syms)

function QEDcore.PhaseSpacePoint(
        p::AbstractProcessDefinition,
        m::AbstractModelDefinition,
        psl::AbstractPhaseSpaceLayout,
        coords::Tuple,
    )

    in_dim = phase_space_dimension(p, m, in_phase_space_layout(psl))
    out_dim = phase_space_dimension(p, m, psl)
    return PhaseSpacePoint(
        p,
        m,
        psl,
        ntuple(i -> coords[i], in_dim),
        ntuple(i -> coords[in_dim + i], out_dim),
    )
end

QEDcore.momenta(psp::AbstractPhaseSpacePoint) = (momenta(psp, Incoming())..., momenta(psp, Outgoing())...)

# PSL for heads-on collision

abstract type AbstractTwoBodyHeadsOnSystem <: AbstractTwoBodyInPhaseSpaceLayout end
QEDbase.phase_space_dimension(proc, model, ::AbstractTwoBodyHeadsOnSystem) = 2 # E, omega

struct PhotonElectronHeadsOnSystem <: AbstractTwoBodyHeadsOnSystem end

Base.broadcastable(psl::PhotonElectronHeadsOnSystem) = Ref(psl)

function QEDbase._build_momenta(
        ::AbstractProcessDefinition,
        ::AbstractPerturbativeModel,
        ::PhotonElectronHeadsOnSystem,
        in_coords::NTuple{2, T}
    ) where {T <: Real}
    E, om = in_coords
    rho = sqrt(E^2 - one(E))

    K = SFourMomentum(om, 0, 0, om)
    P = SFourMomentum(E, 0, 0, -rho)

    return (P, K)
end
