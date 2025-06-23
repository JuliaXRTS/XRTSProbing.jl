module QEDprobing

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


### In-Medium modifications

# temperature
export FiniteTemperature, ZeroTemperature

# matter model
export AbstractMatterModel

# Electron system
export AbstractElectronSystem
export temperature, electron_density, imag_dynamic_response, real_dynamic_response
export fermi_wave_vector,
    fermi_energy, beta, betabar, dynamic_response, dynamic_structure_factor

export AbstractProperElectronSystem
export AbstractInteractingElectronSystem
export proper_electron_system, screening

# screening
export AbstractScreening, Screening, NoScreening
export dielectric_function,
    pseudo_potential, local_field_correction, local_effective_potential
export AbstractPseudoPotential, CoulombPseudoPotential
export AbstractLocalFieldCorrection, NoLocalFieldCorrection

# concrete electron systems
export IdealElectronSystem
export AbstractResponseApproximation, NoApprox, NonDegenerated, Degenerated
export response_approximation
export InteractingElectronSystem

# electron distributions
export EnergyBasedElectronDistribution
export MaxellElectronEnergyDistribution

### Probing Setup
export ProbingSetup
export ElectronProbingSetup
export medium, background_field, kinematic_cuts, differential_cross_section_setup

### Utils

# constants

export HBARC,
    HBARC_eV_ANG,
    ELECTRONMASS,
    ALPHA,
    ALPHA_SQUARE,
    ELEMENTARY_CHARGE_SQUARED,
    ELEMENTARY_CHARGE,
    HARTREE,
    BOHR_RADIUS_ANG

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

using QEDbase
using QEDcore
using QEDprocesses
using QEDevents
using QEDfields

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

include("constants.jl")
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

include("medium/utils.jl")
include("medium/temperature.jl")


include("medium/interface.jl")
include("medium/generic.jl")
include("medium/electron_system/utils.jl")
include("medium/electron_system/interface.jl")
include("medium/electron_system/generic.jl")
include("medium/electron_system/ideal/approximations/interface.jl")
include("medium/electron_system/ideal/approximations/no_approx.jl")
include("medium/electron_system/ideal/approximations/non_degenerated.jl")
include("medium/electron_system/ideal/approximations/degenerated.jl")
include("medium/electron_system/ideal/utils.jl")
include("medium/electron_system/ideal/interface.jl")
include("medium/electron_system/ideal/generic.jl")
include("medium/electron_system/ideal/impl.jl")
include("medium/electron_system/interacting/screening/interface.jl")
include("medium/electron_system/interacting/screening/generic.jl")
include("medium/electron_system/interacting/screening/impl.jl")
include("medium/electron_system/interacting/interface.jl")
include("medium/electron_system/interacting/generic.jl")
include("medium/electron_system/interacting/impl.jl")

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
