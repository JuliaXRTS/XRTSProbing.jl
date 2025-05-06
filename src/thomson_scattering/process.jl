struct Thomson{InElectronSpin, InPhotonPol, OutElectronSpin, OutPhotonPol} <:
    AbstractProcessDefinition where {
        InElectronSpin <: AbstractSpin,
        InPhotonPol <: AbstractPolarization,
        OutElectronSpin <: AbstractSpin,
        OutPhotonPol <: AbstractPolarization,
    }
    in_spin::InElectronSpin
    in_pol::InPhotonPol

    out_spin::OutElectronSpin
    out_pol::OutPhotonPol
end

function Thomson()
    return Thomson(AllSpin(), AllPol(), AllSpin(), AllPol())
end

function Thomson(in_pol::AbstractPolarization)
    return Thomson(AllSpin(), in_pol, AllSpin(), AllPol())
end
function Thomson(in_pol::AbstractPolarization, out_pol::AbstractPolarization)
    return Thomson(AllSpin(), in_pol, AllSpin(), out_pol)
end

function QEDbase.incoming_particles(::Thomson)
    return (Electron(), Photon())
end

function QEDbase.outgoing_particles(::Thomson)
    return (Electron(), Photon())
end

function QEDbase.incoming_spin_pols(proc::Thomson)
    return (proc.in_spin, proc.in_pol)
end

function QEDbase.outgoing_spin_pols(proc::Thomson)
    return (proc.out_spin, proc.out_pol)
end

function _spin_or_pol(process::Thomson, ::Electron, ::Incoming)
    return process.in_spin
end

function _spin_or_pol(process::Thomson, ::Electron, ::Outgoing)
    return process.out_spin
end

function _spin_or_pol(process::Thomson, ::Photon, ::Incoming)
    return process.in_pol
end

function _spin_or_pol(process::Thomson, ::Photon, ::Outgoing)
    return process.out_pol
end

function Base.show(io::IO, ::Thomson)
    print(io, "one-photon Thomson scattering")
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", proc::Thomson)
    println(io, "one-photon Thomson scattering")
    println(io, "    incoming: electron ($(proc.in_spin)), photon ($(proc.in_pol))")
    println(io, "    outgoing: electron ($(proc.out_spin)), photon ($(proc.out_pol))")
    return nothing
end
# TODO: extent this using PSL
function QEDbase.in_phase_space_dimension(proc::Thomson, ::PerturbativeQED)
    return 1 # usually energy of the init photon
end

# TODO: extent this using PSL
function QEDbase.out_phase_space_dimension(proc::Thomson, ::PerturbativeQED)
    return 2 # usually polar and azimuthal angle of the scattered photon
end
