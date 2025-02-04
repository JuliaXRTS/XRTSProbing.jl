module QEDprobing

export PertQED
export FlatPhaseSpaceSampler, randmom

export ThomsonScattering
export ElabPhotonSphSystem

using Random
using QEDbase
using QEDcore
using QEDprocesses
using QEDevents

function hello_world()
    return "Hello, World!"
end

include("patch_QEDprocesses.jl")
include("patch_QEDevents.jl")

include("constants.jl")

include("thomson_scattering/process.jl")
include("thomson_scattering/kinematics.jl")
include("thomson_scattering/cross_section.jl")
include("thomson_scattering/sampler.jl")

end
