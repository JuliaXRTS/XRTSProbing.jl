module XRTSProbing

### Hard scattering

# photon distributions
export XAxis, YAxis, ZAxis, UniformAxis
export GaussianPhotonDist, energy_mean, energy_width, k_vec_axis
export DistributionBasedField, energy_spectrum

# compute models (remove after patch_QEDcore is merged)
export PertQED
export FlatPhaseSpaceSampler, randmom

# scattering processes
export Thomson

# phase space layouts
export PhotonSphericalLayout, Elastic, InElastic, kinematic_mode, is_elastic

# kinematic cuts
export KinematicCuts, is_within_cuts, all_within_cuts
export degree_of_freedom

# electron distributions
export EnergyBasedElectronDistribution
export MaxellElectronEnergyDistribution

### Probing Setup
export ProbingSetup
export ElectronProbingSetup
export medium, background_field, kinematic_cuts, differential_cross_section_setup

### Utils

# compute setups
export AbstractComputationSetup, AbstractProcessSetup
export compute, scattering_process, physical_model

export DifferentialCrossSection
export cached_coords

# remove if patch_QEDcore is resolved
export coordinate_symbol, coordinate_symbols, coord_index
export PhotonElectronHeadsOnSystem

### Event Generation

# proposal generation
export VegasGrid, VegasProposal
export nbins, extent, nodes, spacing
export uniform_vegas_nodes
export train!

# maximum finding
export QuantileReductionMethod

# event generation
export EventGenerator
export generate_event, generate_events

# events container
export Event

using Reexport

using QEDbase
using QEDcore
using QEDprocesses
using QEDevents
using QEDfields
@reexport using ElectronicStructureModels

using Random
using StaticArrays
using LogExpFunctions
using SpecialFunctions
using QuadGK
using TupleTools
using Distributions
using Unitful
using HDF5
using LinearAlgebra

function hello_world()
    return "Hello, World!"
end

include("patch_QEDfields.jl")
include("patch_QEDprocesses.jl")
include("patch_QEDevents.jl")
include("patch_QEDcore.jl")
include("patch_QEDbase.jl")

include("kinematic_cuts.jl")
include("utils.jl")
include("units.jl")

include("fields/utils.jl")
include("fields/photon_dists.jl")
include("fields/types.jl")

include("hard_scattering/thomson_scattering/process.jl")
include("hard_scattering/thomson_scattering/kinematics.jl")
include("hard_scattering/thomson_scattering/cross_section.jl")
include("hard_scattering/thomson_scattering/sampler.jl")

include("electron_dist/types.jl")
include("electron_dist/energy_dist.jl")


include("setups/interface.jl")
include("setups/generic.jl")
include("setups/impl/differential_cross_section.jl")
include("setups/impl/probing.jl")
include("setups/impl/probing_elec.jl")


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
