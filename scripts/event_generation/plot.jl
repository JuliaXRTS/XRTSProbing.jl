using QEDprobing
using QEDcore
using QEDprocesses
using CairoMakie
using PairPlots

# differential cross section
PROC = Thomson()
MODEL = PerturbativeQED()
PSL = PhotonSphericalLayout(PhotonElectronHeadsOnSystem())

event_list = QEDprobing.load("data.h5", PROC, MODEL, PSL)


## projections
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

## plotting

om_sample = energy_transfer.(event_list)
q_sample = three_momentum_transfer.(event_list)


table = (;
    om = om_sample,
    q = q_sample,
)

fig = pairplot(
    table => (PairPlots.Hist(colormap = :jet), PairPlots.MarginDensity(color = :orange)),
    bins = Dict(
        :om => range(minimum(om_sample), maximum(om_sample), 100),
        :q => range(minimum(q_sample), maximum(q_sample), 100),
    )
)

save("en_mom_transfer.pdf", fig)
