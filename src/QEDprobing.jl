module QEDprobing

export PertQED
export FlatPhaseSpaceSampler, randmom

export Thomson
export PhotonSphericalLayout, Elastic, InElastic, kinematic_mode, is_elastic

# structure factor
export temperature, betabar
export LindhardDRF


export DifferentialCrossSection, DifferentialCrossSectionCached
export cached_coords
# remove if patch_QEDcore is resolved
export coordinate_symbol, coordinate_symbols, coord_index

export VegasGrid, VegasProposal
export nbins, extent, nodes, spacing
export uniform_vegas_nodes
export train!

export QuantileReductionMethod

export EventGenerator
export generate_event, generate_events


using Random
using QEDbase
using QEDcore
using QEDprocesses
using QEDevents
using StaticArrays
using LogExpFunctions
using QuadGK
using TupleTools

function hello_world()
    return "Hello, World!"
end

include("patch_QEDprocesses.jl")
include("patch_QEDevents.jl")
include("patch_QEDcore.jl")
include("patch_QEDbase.jl")

include("constants.jl")
include("utils.jl")

include("thomson_scattering/process.jl")
include("thomson_scattering/kinematics.jl")
include("thomson_scattering/cross_section.jl")
include("thomson_scattering/sampler.jl")

include("structure_factor/abstract.jl")
include("structure_factor/dynamic/electron_electron/utils.jl")
include("structure_factor/dynamic/electron_electron/types.jl")
include("structure_factor/dynamic/electron_electron/response_function/utils.jl")
include("structure_factor/dynamic/electron_electron/response_function/lindhard.jl")

include("setup.jl")
include("events.jl")

include("proposal/vegas/utils.jl")
include("proposal/vegas/type.jl")
include("proposal/vegas/access.jl")
include("proposal/vegas/map.jl")
include("proposal/vegas/refine.jl")
include("proposal/vegas/training.jl")
include("proposal/vegas/sampler.jl")

include("max_finder/types.jl")
include("max_finder/findmax.jl")
include("max_finder/naive.jl")
include("max_finder/quantile_reduction.jl")

include("sampler/generate.jl")
include("sampler/filter.jl")
include("sampler/select.jl")
include("sampler/finalize.jl")
include("sampler/types.jl")



end
