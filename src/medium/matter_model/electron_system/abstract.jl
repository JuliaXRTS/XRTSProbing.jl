
### general electron systems

abstract type AbstractElectronSystem <: AbstractMatterModel end

"""

    temperature(::AbstractElectronSystem)

Interface function: return the temperature associated with the electron system.
"""
function temperature end

"""

    electron_density(::AbstractElectronSystem)

Interface function: return the density of electrons associated with the electron system.
"""
function electron_density end

"""

    imag_dynamic_response(::AbstractElectronSystem,om_q::NTuple{2,T})

Interface function: return the value of the imaginary part of the dynamic response function associated with the electron system.
"""
function imag_dynamic_response end

"""

    real_dynamic_response(::AbstractElectronSystem,om_q::NTuple{2,T})

Interface function: return the value of the real part of the dynamic response function associated with the electron system.
"""
function real_dynamic_response end

### proper (non-interacting) electron system
# - electron system without screening

abstract type AbstractProperElectronSystem <: AbstractElectronSystem end


### interacting electron systems
# - electron system with screening

abstract type AbstractInteractingElectronSystem <: AbstractElectronSystem end

"""

    proper_electron_system(::AbstractInteractingElectronSystem)

Interface function: return the proper electron system associated with the interacting system.
"""
function proper_electron_system end

"""

    screening(::AbstractInteractingElectronSystem)

Interface function: return the screening associated with the interacting system.
"""
function screening end

## general screening

abstract type AbstractScreening end
abstract type AbstractPseudoPotential end
abstract type AbstractLocalFieldCorrection end

"""

    pseudo_potential(::AbstractScreening, om_q::Tuple)

Interface function: return value of the pseudo potential for interacting electron systems
"""
function pseudo_potential end

"""

    local_field_correction(::AbstractScreening, om_q::Tuple)

Interface function: return value of the local field correction for interacting electron systems
"""
function local_field_correction end
