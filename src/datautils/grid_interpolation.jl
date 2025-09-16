########
# extrapolation methods
########
abstract type AbstractInterpolationMethod end
struct InterpolExtrapol <: AbstractInterpolationMethod end
struct InterpolEndValue <: AbstractInterpolationMethod end
struct NearestInput <: AbstractInterpolationMethod end
struct BelowInput <: AbstractInterpolationMethod end
struct AboveInput <: AbstractInterpolationMethod end

"""
    GridInterpolant(
        cols::T1,
        rows::T2,
        values::T3
    ) where {T1<:AbstractVector{T},T2<:AbstractVector{T},T3<:AbstractMatrix{T}} where {T<:Real}

TBW
"""
struct GridInterpolant{T1, T2, T3, M}
    rows::T1
    cols::T2
    values::T3
    method::M

    function GridInterpolant(
            cols::T1, rows::T2, values::T3, method::M
        ) where {
            T1 <: AbstractVector{T}, T2 <: AbstractVector{T}, T3 <: AbstractMatrix{T}, M <: AbstractInterpolationMethod,
        } where {T <: Real}
        len_arg_rows = length(rows)
        len_arg_cols = length(cols)
        len_val_rows, len_val_cols = size(values)

        len_arg_cols == len_val_cols || throw(
            ArgumentError(
                "Length of column arguments <$len_arg_cols> does not match the the first dimension of the values matrix <$len_val_cols>.",
            ),
        )

        len_arg_rows == len_val_rows || throw(
            ArgumentError(
                "Length of row arguments <$len_arg_rows> does not match the the second dimension of the values matrix <$len_val_rows>.",
            ),
        )

        # TODO: implement check, if rows and cols are in ascending order

        return new{T1, T2, T3, M}(rows, cols, values, method)
    end

    # TODO: consider implementing a second constructor for external methods, e.g. comming
    # from Interpolations.jl (in this case, method could just store the interpolant from
    # there)
end

@inline function (grid::GridInterpolant)(col, row)
    return _evaluate(grid, grid.method, col, row)
end

@inline interpolate(cols, rows, vals, method) = GridInterpolant(cols, rows, vals, method)

# TODO:
# - change signature of _evaluate to _evaluate(GLT, method, column,row)
# - check behaviour for sampling outside of the _evaluate
# - add some I/O for _evaluate tables (simple saving and loading)


"""
    _evaluate(
    method::InterpolExtrapol,
    column::Real,
    row::Real,
    table::GridInterpolant)

TBW
"""
function _evaluate(grid::GridInterpolant, method::InterpolExtrapol, x::Real, y::Real)
    # Find the lower indices using binary search
    row_lo = searchsortedfirst(grid.rows, y)
    col_lo = searchsortedfirst(grid.cols, x)

    # Adjust indices to avoid going out of bounds
    if row_lo == length(grid.rows)
        row_lo -= 1
    end
    if col_lo == length(grid.cols)
        col_lo -= 1
    end

    row_hi = row_lo + 1
    col_hi = col_lo + 1

    # Compute interpolation weights (normalized position within cell)
    t_row = (y - grid.rows[row_lo]) / (grid.rows[row_hi] - grid.rows[row_lo])
    t_col = (x - grid.cols[col_lo]) / (grid.cols[col_hi] - grid.cols[col_lo])

    # Values at the four corners of the cell
    v_ll = grid.values[row_lo, col_lo]  # lower-left
    v_lh = grid.values[row_lo, col_hi]  # lower-right
    v_hl = grid.values[row_hi, col_lo]  # upper-left
    v_hh = grid.values[row_hi, col_hi]  # upper-right

    # Bilinear interpolation
    v_bottom = v_ll + t_col * (v_lh - v_ll)
    v_top = v_hl + t_col * (v_hh - v_hl)
    return v_bottom + t_row * (v_top - v_bottom)
end

"""
    _evaluate(
    method::InterpolEndValue,
    column::Real,
    row::Real,
    table::GridInterpolant)

TBW
"""
function _evaluate(grid::GridInterpolant, method::InterpolEndValue, x::Real, y::Real)
    rows = grid.rows
    cols = grid.cols
    vals = grid.values

    # --- Find lower cell indices ---
    row_lo = searchsortedfirst(rows, y)
    col_lo = searchsortedfirst(cols, x)

    # Clamp lower indices so row_hi, col_hi are in bounds
    row_lo = min(max(row_lo, 1), length(rows) - 1)
    col_lo = min(max(col_lo, 1), length(cols) - 1)

    row_hi = row_lo + 1
    col_hi = col_lo + 1

    # --- Handle column interpolation factor ---
    if x >= cols[col_hi]
        # Outside on the right: snap to right column
        t_col = 0.0
        col_lo = col_hi
        col_delta_low = 0.0
    elseif x <= cols[col_lo]
        # Outside on the left: snap to left column
        t_col = 0.0
        col_delta_low = 0.0
    else
        # Inside cell horizontally: do normal interpolation
        t_col = (x - cols[col_lo]) / (cols[col_hi] - cols[col_lo])
        col_delta_low = vals[row_lo, col_hi] - vals[row_lo, col_lo]
    end

    # --- Handle row interpolation factor ---
    if y >= rows[row_hi]
        # Outside above: snap to upper row
        t_row = 0.0
        row_lo = row_hi
        row_delta = 0.0
    elseif y <= rows[row_lo]
        # Outside below: snap to lower row
        t_row = 0.0
        row_delta = 0.0
    else
        # Inside cell vertically: compute normal row interpolation
        col_delta_high = vals[row_hi, col_hi] - vals[row_hi, col_lo]
        t_row = (y - rows[row_lo]) / (rows[row_hi] - rows[row_lo])
        row_delta = (vals[row_hi, col_lo] + col_delta_high * t_col) -
            (vals[row_lo, col_lo] + col_delta_low * t_col)
    end

    return vals[row_lo, col_lo] + t_col * col_delta_low + t_row * row_delta
end


"""
    _evaluate(
    method::NearestInput,
    column::Real,
    row::Real,
    table::GridInterpolant)

TBW
"""
function _evaluate(grid::GridInterpolant, method::NearestInput, x::Real, y::Real)
    rows = grid.rows
    cols = grid.cols
    vals = grid.values

    # --- Find lower cell indices (binary search) ---
    row_lo = searchsortedfirst(rows, y)
    col_lo = searchsortedfirst(cols, x)

    # Clamp indices to valid range
    row_lo = min(max(row_lo, 1), length(rows) - 1)
    col_lo = min(max(col_lo, 1), length(cols) - 1)

    row_hi = row_lo + 1
    col_hi = col_lo + 1

    # --- Determine nearest neighbor row ---
    nearest_row =
    if y < (rows[row_lo] + rows[row_hi]) / 2
        row_lo
    else
        row_hi
    end

    # --- Determine nearest neighbor column ---
    nearest_col =
    if x < (cols[col_lo] + cols[col_hi]) / 2
        col_lo
    else
        col_hi
    end

    return vals[nearest_row, nearest_col]
end

"""
    _evaluate(
    method::BelowInput,
    column::Real,
    row::Real,
    table::GridInterpolant)

TBW
"""
function _evaluate(grid::GridInterpolant, method::BelowInput, x::Real, y::Real)
    rows = grid.rows
    cols = grid.cols
    vals = grid.values

    # --- Find lower cell indices (binary search) ---
    row_lo = searchsortedfirst(rows, y)
    col_lo = searchsortedfirst(cols, x)

    # Clamp indices so they never exceed grid bounds
    row_lo = min(max(row_lo, 1), length(rows))
    col_lo = min(max(col_lo, 1), length(cols))

    return vals[row_lo, col_lo]
end


"""
    _evaluate(
    method::AboveInput,
    column::Real,
    row::Real,
    table::GridInterpolant)

TBW
"""
function _evaluate(grid::GridInterpolant, method::AboveInput, x::Real, y::Real)
    rows = grid.rows
    cols = grid.cols
    vals = grid.values

    # --- Find lower cell indices (binary search) ---
    row_lo = searchsortedfirst(rows, y)
    col_lo = searchsortedfirst(cols, x)

    # Choose upper index if not exactly on a grid node
    row_hi = (rows[row_lo] == y) ? row_lo : min(row_lo + 1, length(rows))
    col_hi = (cols[col_lo] == x) ? col_lo : min(col_lo + 1, length(cols))

    return vals[row_hi, col_hi]
end
