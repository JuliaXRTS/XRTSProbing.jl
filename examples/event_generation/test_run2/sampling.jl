include("init.jl")

DATADIR = "data"


@info "generating $N events with batchsize $BATCH_SIZE ..."
event_list = generate_events(RNG, EG, N, BATCH_SIZE)
@info "done."

filename = "events.h5"
save_path = joinpath(DATADIR, filename)
QEDprobing.save(event_list, save_path)
@info "Saved event list in $save_path"

filename = "photon_en_mag.h5"
save_path = joinpath(DATADIR, filename)
QEDprobing.save_photon_en_mag(event_list, save_path, keV = true)
@info "Saved photon energy and magnitude in $save_path"
