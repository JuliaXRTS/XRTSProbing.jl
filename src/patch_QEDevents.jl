# quick and dirty implementation of a flat phase space generator (based on RAMBO)

struct FlatPhaseSpaceSampler{
    P<:AbstractProcessDefinition,
    M<:AbstractModelDefinition,
    PSL<:FlatPhaseSpaceLayout,
    CM<:QEDcore.AbstractCoordinateMap,
}
    proc::P
    model::M
    psl::PSL
    cmap::CM

    # currently we only support cached coordinate maps
    function FlatPhaseSpaceSampler(
        proc::P,
        model::M,
        in_psl::PSL,
        in_coords::NTuple{N,T},
    ) where {
        P<:AbstractProcessDefinition,
        M<:AbstractModelDefinition,
        PSL<:AbstractInPhaseSpaceLayout,
        N,
        T<:Real,
    }
        psl = FlatPhaseSpaceLayout(in_psl)
        cmap = CoordinateMapCached(proc, model, psl, in_coords)
        return new{P,M,typeof(psl),typeof(cmap)}(proc, model, psl, cmap)
    end
end

function randmom(
    rng::AbstractRNG,
    smpl::FlatPhaseSpaceSampler{P,M,PSL,CM},
) where {P,M,PSL,CM<:CoordinateMapCached}
    in_moms = smpl.cmap.in_moms
    nin = number_incoming_particles(smpl.proc)
    nout = number_outgoing_particles(smpl.proc)
    out_rnd_coords = NTuple{4 * nout,Float64}(rand(rng, 4 * nout))
    out_moms = NTuple{nout,SFourMomentum}(smpl.cmap(out_rnd_coords))

    #return NTuple{nin + nout,SFourMomentum}((in_moms..., out_moms...))
    return in_moms, out_moms
end

randmom(smpl::FlatPhaseSpaceSampler) = randmom(Random.default_rng(), smpl)
randmom(smpl::FlatPhaseSpaceSampler, n::Int) = randmom(Random.default_rng(), smpl, n)

function randmom(rng::AbstractRNG, smpl::FlatPhaseSpaceSampler, n::Int)
    nin = number_incoming_particles(smpl.proc)
    nout = number_outgoing_particles(smpl.proc)
    res = Vector{NTuple{nout + nin,SFourMomentum}}(undef, n)

    @inbounds for i in eachindex(res)
        res[i] = randmom(rng, smpl)
    end

    return res
end
