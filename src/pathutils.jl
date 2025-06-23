# utility functions for the data
#
# wording:
# - 'datadir' refers to the path within the data folder
# - `readdir` refers to the content of the folder behind a given path
# - `datapaths` refers to a collection of paths to data files and folders

datadir_raw(subdir...) = datadir("raw", subdir...)
datadir_raw_ext(subdir...) = datadir("raw", "external", subdir...)
datadir_processed(subdir...) = datadir("processed", subdir...)
datadir_processed_ext(subdir...) = datadir("processed", "external", subdir...)

readdir_raw_data(subdir...) = readdir(datadir_raw(subdir...))
readdir_raw_ext_data(subdir...) = readdir(datadir_raw_ext(subdir...))
readdir_processed_data(subdir...) = readdir(datadir_processed(subdir...))
readdir_processed_ext_data(subdir...) = readdir(datadir_processed_ext(subdir...))

"""
Exchanges the file extension of `file_name` with `new_ext`. No additional checks are performed.
"""
function _file_ext_to(file_name, new_ext)
    tmp = split(file_name, ".")
    tmp[end] = new_ext
    return join(tmp, ".")
end
