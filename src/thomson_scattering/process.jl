# TODO: introduce Incident photon polarization: PolX, PolY, AllPol (which is the current
# case)
struct ThomsonScattering <: AbstractProcessDefinition end

function QEDprocesses.incoming_particles(::ThomsonScattering)
    return (Electron(), Photon())
end

function QEDprocesses.outgoing_particles(::ThomsonScattering)
    return (Electron(), Photon())
end

# TODO: extent this using PSL
function QEDbase.in_phase_space_dimension(proc::ThomsonScattering, ::PerturbativeQED)
    return 1 # usually energy of the init photon
end

# TODO: extent this using PSL
function QEDbase.out_phase_space_dimension(proc::ThomsonScattering, ::PerturbativeQED)
    return 2 # usually polar and azimuthal angle of the scattered photon
end
