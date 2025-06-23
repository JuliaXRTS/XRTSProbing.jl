# I/O for raw MCSS files

"""
    _read_MCSS_q_array_file(q_array_path)

Reads and parsed the q-array file at `q_array_path` and returns a dataframe with the respective theta and q values.
"""
function _read_MCSS_q_array_file(q_array_path)
    stat_file_names = String[]
    out_file_names = String[]
    thetas = Float64[]
    theta_units = String[]
    wavenumbers = Float64[]
    wavenumber_units = String[]

    open(q_array_path) do file
        for ln in eachline(file)
            sp = split(ln, " ")
            stat_file_name = chop(sp[1]; tail = 1)
            push!(stat_file_names, stat_file_name)
            push!(out_file_names, _file_ext_to(stat_file_name, "csv"))
            push!(thetas, parse(Float64, sp[2]))
            push!(theta_units, chop(sp[3]; tail = 1))
            push!(wavenumbers, parse(Float64, sp[8]))
            push!(wavenumber_units, sp[9])
        end
    end
    return DataFrame(
        :stat_file_name => stat_file_names,
        :out_file_name => out_file_names,
        :theta => thetas,
        :theta_unit => theta_units,
        :wavenumber => wavenumbers,
        :wavenumber_unit => wavenumber_units,
    )
end

"""
    _read_MCSS_stat_file(stat_path, stat_file_name)

Read and parse a MCSS stat file with the filename `stat_file_name` at the absolute path `stat_path`.
"""
function _read_MCSS_stat_file(stat_path, stat_file_name)
    stat_file_path = joinpath(stat_path, stat_file_name)
    file_content = String[]

    # reading the file content
    open(stat_file_path) do file
        for ln in eachline(file)
            push!(file_content, ln)
        end
    end

    # regex for calculated values: "The <key> is: <value>"
    r = r"The ([A-Z]|[a-z]| ){1,} is:"

    # iterate the lines
    line_pointer = 1
    init_param_dict = Dict{String, Dict}()
    calc_val_dict = Dict{String, String}()
    while line_pointer <= length(file_content)
        line = file_content[line_pointer]

        # find and parse init parameter blocks: "$<block name> <block> $END"
        # where <block> is a list which matches: "<key> = <value>"
        if occursin("\$", line) && !occursin("END", line)
            block_name = chop(line; head = 1, tail = 0)
            block_pointer = line_pointer + 1
            block_run = true
            #println("block start")
            #println("----- ",block_name)
            block_dict = Dict{String, String}()
            while block_run
                block_line = file_content[block_pointer]
                if occursin("\$END", block_line)
                    block_run = false
                    #println("block end")
                else
                    key_val = Pair(strip.(split(block_line, " = "))...)
                    #println(key_val)
                    push!(block_dict, key_val)
                    block_pointer += 1
                end
            end
            push!(init_param_dict, block_name => block_dict)
            line_pointer = block_pointer
        end

        # find and parse calculated values: "The <key> is: <value>"
        if !isnothing(match(r, line))
            #println(line)
            pre, val = strip.(split(line, ": "))
            key = join(split(pre, " ")[2:(end - 1)], "_")
            push!(calc_val_dict, Pair(key, val))
        end
        line_pointer += 1
    end
    return Dict("INIT_PARAMS" => init_param_dict, "CALC_VALS" => calc_val_dict)
end

"""
    _collect_MCSS_stats(dpaths::ExtDataPaths, df_qarray)

Collect MCSS stats data from the stat-file names in `df_qarray` located at `dpaths`. Return dict with all stats and calculated data as dataframes.
"""
function _collect_MCSS_stats(dpaths; df_qarray = _read_MCSS_q_array_file(dpaths.qarray))
    df_calc_data = DataFrame()
    df_DETECTOR = DataFrame()
    df_MODELS = DataFrame()
    df_FILE = DataFrame()
    df_PLASMA = DataFrame()
    df_SETUP = DataFrame()
    df_PROBE = DataFrame()

    for row in eachrow(df_qarray)
        #println(row)
        th = row.theta
        stat_file_name = row.stat_file_name
        single_stat_dict = _read_MCSS_stat_file(dpaths.stat, stat_file_name)

        calc_val_dict = copy(single_stat_dict["CALC_VALS"])
        push!(calc_val_dict, "theta" => "$th")
        append!(df_calc_data, calc_val_dict)

        init_param_dict = copy(single_stat_dict["INIT_PARAMS"])

        DETECTOR_dict = copy(init_param_dict["DETECTOR"])
        push!(DETECTOR_dict, "theta" => "$th")
        append!(df_DETECTOR, DETECTOR_dict)

        MODELS_dict = copy(init_param_dict["MODELS"])
        push!(MODELS_dict, "theta" => "$th")
        append!(df_MODELS, MODELS_dict)

        FILE_dict = copy(init_param_dict["FILE"])
        push!(FILE_dict, "theta" => "$th")
        append!(df_FILE, FILE_dict)

        PLASMA_dict = copy(init_param_dict["PLASMA"])
        push!(PLASMA_dict, "theta" => "$th")
        append!(df_PLASMA, PLASMA_dict)

        SETUP_dict = copy(init_param_dict["SETUP"])
        push!(SETUP_dict, "theta" => "$th")
        append!(df_SETUP, SETUP_dict)

        PROBE_dict = copy(init_param_dict["PROBE"])
        push!(PROBE_dict, "theta" => "$th")
        append!(df_PROBE, PROBE_dict)
    end

    df_dict = Dict(
        "CALC_DATA" => mapcols(col -> parse_entry.(col), df_calc_data),
        "DETECTOR" => mapcols(col -> parse_entry.(col), df_DETECTOR),
        "MODELS" => mapcols(col -> parse_entry.(col), df_MODELS),
        "FILE" => mapcols(col -> parse_entry.(col), df_FILE),
        "PLASMA" => mapcols(col -> parse_entry.(col), df_PLASMA),
        "SETUP" => mapcols(col -> parse_entry.(col), df_SETUP),
        "PROBE" => mapcols(col -> parse_entry.(col), df_PROBE),
    )

    return df_dict
end

"""
    _read_MCSS_out_file(path, out_file_name)

Read MCSS outfile `out_file_name` located in `out_path` and parse its content into a dataframe.
"""
function _read_MCSS_out_file(out_path, out_file_name)
    out_file_path = joinpath(out_path, out_file_name)
    header = open(out_file_path) do file
        first(eachline(file))
    end
    sp_header = strip.(split(header[2:end], "]"))[1:(end - 1)] .* "]"
    df = CSV.read(out_file_path, DataFrame; header = 0, skipto = 2)
    rename!(df, sp_header)
    return df
end

"""
    _collect_MCSS_out(dpaths::ExtDataPaths)

Collect values from the MCSS out-files located at `dpaths`. Return dict with `q`, `omega`, `Ps`, `Ps_el`, `Ps_bf`, `Ps_ff`. The omega matrix will be projected on a single vector if possible, i.e. if all columns are the same.
"""
function _collect_MCSS_out(dpaths)
    df_qarray = _read_MCSS_q_array_file(dpaths.qarray)

    theta = df_qarray[!, :theta]
    theta_unit = df_qarray[!, :theta_unit][1]
    wavenumber = df_qarray[!, :wavenumber]
    wavenumber_unit = df_qarray[!, :wavenumber_unit][1]

    len_q = size(df_qarray, 1)
    file_name = df_qarray[!, :out_file_name][1]
    data = QEDprobing._read_MCSS_out_file(dpaths.results, file_name)
    omega_1 = data[!, "E [eV]"]
    omega_1_len = size(omega_1, 1)

    omega_mat = Matrix{Float64}(undef, omega_1_len, len_q)
    Ps_mat = Matrix{Float64}(undef, omega_1_len, len_q)
    Ps_el_mat = Matrix{Float64}(undef, omega_1_len, len_q)
    Ps_bf_mat = Matrix{Float64}(undef, omega_1_len, len_q)
    Ps_ff_mat = Matrix{Float64}(undef, omega_1_len, len_q)

    for (i, row) in enumerate(eachrow(df_qarray))
        file_name = row[:out_file_name]
        data = QEDprobing._read_MCSS_out_file(dpaths.results, file_name)
        omega_mat[:, i] .= reverse(data[!, "E [eV]"])

        Ps_mat[:, i] .= reverse(data[!, "P_{s} [Arb.]"])
        Ps_el_mat[:, i] .= reverse(data[!, "P_{s}^{el} [Arb.]"])
        Ps_bf_mat[:, i] .= reverse(data[!, "P_{s}^{bf} [Arb.]"])
        Ps_ff_mat[:, i] .= reverse(data[!, "P_{s}^{ff} [Arb.]"])
    end

    omega_out = _project_to_1d(omega_mat)

    # TODO: consider using AxisArrays for this!
    return Dict(
        "q" => wavenumber,
        "omega" => omega_out,
        "Ps" => Ps_mat,
        "Ps_el" => Ps_el_mat,
        "Ps_bf" => Ps_bf_mat,
        "Ps_ff" => Ps_ff_mat,
    )
end
