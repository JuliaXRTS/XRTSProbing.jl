# path utility specific for MCSS data
#

datadir_raw_MCSS(subdir...) = datadir_raw_ext("MCSS", subdir...)
readdir_raw_MCSS(subdir...) = readdir(datadir_raw_MCSS(subdir...))

datadir_processed_MCSS(subdir...) = datadir_processed_ext("MCSS", subdir...)
readdir_processed_MCSS(subdir...) = readdir_processed_data("MCSS", subdir...)

datadir_raw_testdata_MCSS(subdir...) = datadir_raw_MCSS("test_data", subdir...)
readdir_raw_testdata_MCSS(subdir...) = readdir(datadir_raw_testdata_MCSS(subdir...))

datadir_processed_testdata_MCSS(subdir...) = datadir_processed_MCSS("test_data", subdir...)
function readdir_processed_testdata_MCSS(subdir...)
    return readdir(datadir_processed_testdata_MCSS(subdir...))
end

# TODO: consider renaming to MCSSDataPaths
Base.@kwdef struct ExtDataPaths
    data
    qarray
    stat
    results
    subdir
    success
end

function Base.show(io::IO, dpath::ExtDataPaths)
    println(io, "ExtDataPaths")
    for field in (:data, :qarray, :stat, :results, :subdir, :success)
        println(io, "\t", string(field), ": ", getfield(dpath, field))
    end
    return
end

"""
Return a dictionary with keys `(:data, :qarray, :stat, :out)` whose values are the respective files/folder.
Throws an error if the data path does not exist and shows warnings if one of the other files/folders do not exist.
"""
function datapaths_raw_MCSS(subdir...; q_array_file_name = "q_array.txt")
    data_path = datadir_raw_MCSS(subdir...)
    q_array_file_path = joinpath(data_path, q_array_file_name)
    stat_path = joinpath(data_path, "stat")
    out_path = joinpath(data_path, "out")

    success = true

    isdir(data_path) || throw(ArgumentError("There is no folder <$data_path>."))
    isfile(q_array_file_path) || begin
        success = false
        @warn(
            "There is no file <$q_array_file_name> in <$data_path>.\n Maybe <$data_path> is not the right data path."
        )
    end
    isdir(stat_path) || begin
        success = false
        @warn(
            "There is no folder <$stat_path>.\n Maybe <$data_path> is not the right data path, or has not the (out, stat) format."
        )
    end
    isdir(out_path) || begin
        sucess = false
        @warn(
            "There is no folder <$out_path>.\n Maybe <$data_path> is not the right data path, or has not the (out, stat) format."
        )
    end
    return ExtDataPaths(;
        data = data_path,
        qarray = q_array_file_path,
        stat = stat_path,
        results = out_path,
        success = success,
        subdir = subdir,
    )
end
