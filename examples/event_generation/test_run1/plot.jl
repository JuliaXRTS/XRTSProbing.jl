using QEDprobing
using QEDcore
using QEDprocesses
using Unitful
using BenchmarkTools
using Random
#using Plots
using StatsPlots
using LaTeXStrings

DATADIR = "data"
PLOTDIR = "plots"

# differential cross section
PROC = Thomson()
MODEL = PerturbativeQED()
PSL = PhotonSphericalLayout(PhotonElectronHeadsOnSystem(XAxis()))

event_list = QEDprobing.load(joinpath(DATADIR, "events.h5"), PROC, MODEL, PSL)

momentum_transfer(ev::Event) = momentum_transfer(ev.psp)
function momentum_transfer(psp::AbstractPhaseSpacePoint)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    out_photon_mom = momentum(psp, Outgoing(), Photon())
    return out_photon_mom - in_photon_mom
end

energy_transfer(ev::Event) = energy_transfer(ev.psp)
function energy_transfer(psp::AbstractPhaseSpacePoint)
    return getE(momentum_transfer(psp))
end

three_momentum_transfer(ev::Event) = three_momentum_transfer(ev.psp)
function three_momentum_transfer(psp::AbstractPhaseSpacePoint)
    return getMag(momentum_transfer(psp))
end

in_photon_energy(ev::Event) = in_photon_energy(ev.psp)
function in_photon_energy(psp::AbstractPhaseSpacePoint)
    in_photon_mom = momentum(psp, Incoming(), Photon())
    return getE(in_photon_mom)
end

out_photon_cth(ev::Event) = out_photon_cth(ev.psp)
function out_photon_cth(psp::AbstractPhaseSpacePoint)
    out_photon_mom = momentum(psp, Outgoing(), Photon())
    return getCosTheta(out_photon_mom)
end

om_sample = energy_transfer.(event_list)
q_sample = three_momentum_transfer.(event_list)
omx_sample = in_photon_energy.(event_list)
cth_sample = out_photon_cth.(event_list)

#P = histogram(q_sample,xlim = (0.0,0.01))
P = histogram(q_sample)
filename = "momentum_transfer.pdf"
save_path = joinpath(PLOTDIR, filename)
savefig(P, save_path)

#P = histogram(om_sample,xlim = (-1e-5,0.0))
P = histogram(om_sample)
filename = "energy_transfer.pdf"
save_path = joinpath(PLOTDIR, filename)
savefig(P, save_path)

P = histogram(cth_sample)
filename = "scattering_angle.pdf"
save_path = joinpath(PLOTDIR, filename)
savefig(P, save_path)


#P = histogram2d(om_sample, q_sample,bins=(40,40),cmap=:plasma,show_empty=true,dpi=600,normalize=:pdf,size=(800,800))
P = marginalhist(
    om_sample,
    q_sample,
    nbins = 100,
    xlab = L"\omega \quad [m_e]",
    ylab = L"q \quad [m_e]",
    cmap = :binary,
    show_empty = true
)
filename = "energy_momentum_transfer.png"
save_path = joinpath(PLOTDIR, filename)
savefig(P, save_path)

P = histogram(omx_sample, nbins = 100, cmap = :magma)
filename = "init_photon_energy.pdf"
save_path = joinpath(PLOTDIR, filename)
savefig(P, save_path)
