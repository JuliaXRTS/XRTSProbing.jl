# I/O for processed MCSS files

"""
    _save_MCSS_out_data(data_dict, dpaths::ExtDataPaths; file_name="data.jld2")

Saves the out data-dict `data_dict` (e.g. return of `_collect_MCSS_out`) into `data/processed/external/MCSS/subdir...`, where `subdir` is given by `dpaths`.

!!! note "safe save"

    If the file to be saved already exists, the existing file will be alternated with a version number and the new data will be saved with the given name.

"""
function _save_MCSS_out_data(data_dict, dpaths; file_name = "data.jld2")
    data_dir = datadir_processed_MCSS(dpaths.subdir...)
    mkpath(data_dir)
    file_path = joinpath(data_dir, file_name)
    return safesave(file_path, data_dict)
end

"""
    _load_MCSS_out_data(dpaths, file_name="data.jld2")

Load the processed MCSS out-file from `data/processed/external/MCSS/subdir...`, where `subdir` is given by `dpaths`. Return a data-dict with q-array, omega-array, and all saved quantities as matrices.
"""
function _load_MCSS_out_data(dpaths; file_name = "data.jld2")
    data_dir = datadir_processed_MCSS(dpaths.subdir...)
    file_path = joinpath(data_dir, file_name)
    return load(file_path)
end
