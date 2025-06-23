"""
    _project_to_1d(arr::AbstractMatrix)

Return first column of a matrix, if all columns are the same. Return the matrix itself otherwise.
"""
function _project_to_1d(arr::AbstractMatrix)
    projectable = true
    col1 = arr[:, 1]
    for col in eachcol(arr)
        if !isapprox(col, col1)
            projectable = false
            break
        end
    end
    if projectable
        return col1
    else
        return arr
    end
end
