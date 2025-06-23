function parse_entry(st)
    #@show st

    m_false = match(r"\.FALSE\.", st)
    if !isnothing(m_false)
        # parse m_false
        return false
    end

    m_true = match(r"\.TRUE\.", st)
    if !isnothing(m_true)
        return true
    end

    m_str = match(r"\".{1,}\"", st)
    if !isnothing(m_str)
        # parse m_str
        return String(chop(st; head = 1, tail = 1))
    end
    #m_number = match(r"(\d{1,}\.\d{1,}E(\+|\-){1}\d{1,})",st)
    #m_number = match(r"\d",st)

    m_number = match(r"^(?![a-zA-z])(\d|-)", st)
    if !isnothing(m_number)
        sp_st = split(st, " ")
        if length(sp_st) == 1
            return _parse_number(st)
        else
            return (_parse_number(sp_st[1]), sp_st[2]) # maybe use Unitful?
        end
    end

    m_array = match(r"^\(.+\)$", st)
    if !isnothing(m_array)
        arr_st = strip.(split(chop(st; head = 1, tail = 1), ","))
        return _parse_number.(arr_st)
    end

    m_element = match(r"([a-zA-Z]+?\(.+\))", st)
    if !isnothing(m_element)
        return String(st)
    end

    throw(ArgumentError("<$st> does not match any default pattern"))
end

function _parse_number(st)
    if occursin(".", st) || occursin("e", st) || occursin("E", st)
        return parse(Float64, st)
    else
        return parse(Int, st)
    end
end
