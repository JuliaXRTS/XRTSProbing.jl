########
# Extrapolation / Interpolation Methods
########

abstract type AbstractInterpolationMethod end

"""
    InterpolExtrapol

Bilinear interpolation allowing linear extrapolation outside the grid.
"""
struct InterpolExtrapol <: AbstractInterpolationMethod end

"""
    InterpolEndValue

Bilinear interpolation clamped to the edge values (no extrapolation).
"""
struct InterpolEndValue <: AbstractInterpolationMethod end

"""
    NearestInput

Nearest-neighbor interpolation: returns the value of the nearest grid node.
"""
struct NearestInput <: AbstractInterpolationMethod end

"""
    BelowInput

Step-wise interpolation returning the value from the grid cell below and to the left.
"""
struct BelowInput <: AbstractInterpolationMethod end

"""
    AboveInput

Step-wise interpolation returning the value from the grid cell above and to the right.
"""
struct AboveInput <: AbstractInterpolationMethod end

########
# GridInterpolant
########

"""
    GridInterpolant(cols, rows, values, method)

Represents a 2D grid interpolant.

# Arguments
- `cols::AbstractVector{T}`: Grid points along the column (x-axis).
- `rows::AbstractVector{T}`: Grid points along the row (y-axis).
- `values::AbstractMatrix{T}`: Values at grid nodes, dimension `(length(rows), length(cols))`.
- `method::AbstractInterpolationMethod`: Interpolation method.
"""
struct GridInterpolant{T1, T2, T3, M}
    rows::T1
    cols::T2
    values::T3
    method::M

    function GridInterpolant(cols::T1, rows::T2, values::T3, method::M) where {
            T <: Real, T1 <: AbstractVector{T}, T2 <: AbstractVector{T}, T3 <: AbstractMatrix{T}, M <: AbstractInterpolationMethod,
        }
        length(rows) == size(values, 1) ||
            throw(ArgumentError("Length of rows does not match first dimension of values"))
        length(cols) == size(values, 2) ||
            throw(ArgumentError("Length of columns does not match second dimension of values"))

        return new{T1, T2, T3, M}(rows, cols, values, method)
    end
end

# Make GridInterpolant callable
@inline function (grid::GridInterpolant)(x, y)
    return _evaluate(grid, grid.method, x, y)
end

@inline interpolate(cols, rows, vals, method) = GridInterpolant(cols, rows, vals, method)

########
# Reusable helper: find cell indices and weights
########
"""
    _find_cell(grid, x, y)

Finds the lower and upper indices in rows and columns and
computes normalized weights within the cell.

Returns `(row_lo, row_hi, t_row, col_lo, col_hi, t_col)`.
"""
function _find_cell(grid::GridInterpolant, x::Real, y::Real)
    rows, cols = grid.rows, grid.cols

    # Lower indices
    row_lo = searchsortedfirst(rows, y)
    col_lo = searchsortedfirst(cols, x)

    # Clamp to ensure row_hi and col_hi exist
    row_lo = min(max(row_lo, 1), length(rows) - 1)
    col_lo = min(max(col_lo, 1), length(cols) - 1)

    row_hi = row_lo + 1
    col_hi = col_lo + 1

    # Normalized weights in [0,1]
    t_row = (y - rows[row_lo]) / (rows[row_hi] - rows[row_lo])
    t_col = (x - cols[col_lo]) / (cols[col_hi] - cols[col_lo])

    return row_lo, row_hi, t_row, col_lo, col_hi, t_col
end

########
# Evaluate implementations
########

# --- Bilinear interpolation with linear extrapolation ---
function _evaluate(grid::GridInterpolant, ::InterpolExtrapol, x::Real, y::Real)
    row_lo, row_hi, t_row, col_lo, col_hi, t_col = _find_cell(grid, x, y)
    v = grid.values

    # Bilinear interpolation
    v_bottom = v[row_lo, col_lo] + t_col * (v[row_lo, col_hi] - v[row_lo, col_lo])
    v_top = v[row_hi, col_lo] + t_col * (v[row_hi, col_hi] - v[row_hi, col_lo])
    return v_bottom + t_row * (v_top - v_bottom)
end

# --- Bilinear interpolation clamped to edge ---
function _evaluate(grid::GridInterpolant, ::InterpolEndValue, x::Real, y::Real)
    rows, cols, vals = grid.rows, grid.cols, grid.values

    # --- Find lower cell indices ---
    row_lo = searchsortedfirst(rows, y)
    col_lo = searchsortedfirst(cols, x)

    # Clamp indices to ensure row_hi and col_hi exist
    row_lo = min(max(row_lo, 1), length(rows) - 1)
    col_lo = min(max(col_lo, 1), length(cols) - 1)
    row_hi = row_lo + 1
    col_hi = col_lo + 1

    # --- Compute horizontal (column) factor ---
    if x >= cols[col_hi]
        t_col = 0.0
        col_lo = col_hi
        col_delta_low = 0.0
    elseif x <= cols[col_lo]
        t_col = 0.0
        col_delta_low = 0.0
    else
        t_col = (x - cols[col_lo]) / (cols[col_hi] - cols[col_lo])
        col_delta_low = vals[row_lo, col_hi] - vals[row_lo, col_lo]
    end

    # --- Compute vertical (row) factor ---
    if y >= rows[row_hi]
        t_row = 0.0
        row_lo = row_hi
        row_delta = 0.0
    elseif y <= rows[row_lo]
        t_row = 0.0
        row_delta = 0.0
    else
        col_delta_high = vals[row_hi, col_hi] - vals[row_hi, col_lo]
        t_row = (y - rows[row_lo]) / (rows[row_hi] - rows[row_lo])
        row_delta = (vals[row_hi, col_lo] + col_delta_high * t_col) -
            (vals[row_lo, col_lo] + col_delta_low * t_col)
    end

    # --- Return clamped bilinear interpolation ---
    return vals[row_lo, col_lo] + col_delta_low * t_col + row_delta * t_row
end

# --- Nearest neighbor ---
function _evaluate(grid::GridInterpolant, ::NearestInput, x::Real, y::Real)
    row_lo, row_hi, _, col_lo, col_hi, _ = _find_cell(grid, x, y)
    rows, cols = grid.rows, grid.cols

    nearest_row = (y < (rows[row_lo] + rows[row_hi]) / 2) ? row_lo : row_hi
    nearest_col = (x < (cols[col_lo] + cols[col_hi]) / 2) ? col_lo : col_hi
    return grid.values[nearest_row, nearest_col]
end

# --- Stepwise below ---
function _evaluate(grid::GridInterpolant, ::BelowInput, x::Real, y::Real)
    row_lo = searchsortedfirst(grid.rows, y)
    col_lo = searchsortedfirst(grid.cols, x)

    # Clamp
    row_lo = min(max(row_lo, 1), length(grid.rows))
    col_lo = min(max(col_lo, 1), length(grid.cols))
    return grid.values[row_lo, col_lo]
end

# --- Stepwise above ---
function _evaluate(grid::GridInterpolant, ::AboveInput, x::Real, y::Real)
    row_lo = searchsortedfirst(grid.rows, y)
    col_lo = searchsortedfirst(grid.cols, x)
    rows, cols = grid.rows, grid.cols

    row_hi = (rows[row_lo] == y) ? row_lo : min(row_lo + 1, length(rows))
    col_hi = (cols[col_lo] == x) ? col_lo : min(col_lo + 1, length(cols))
    return grid.values[row_hi, col_hi]
end
