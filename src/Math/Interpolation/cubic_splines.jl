export InterpCubicSplines

"""
    InterpCubicSplines{T, N} <: AbstractInterpolationMethod

Type storing a cubic spline nodes and coefficients. `T` is the spline data type and 
`N` is the spline dimension (i.e., the number of interpolated functions). 

### Fields 
- `n` -- Number of node points.
- `xn` -- Interpolated node points.
- `yn` -- Node points function values. 
- `c` --  Spline polynomials coefficients. 
- `type` -- Boundary conditions type.
"""
struct InterpCubicSplines{T,N} <: jMath.AbstractInterpolationMethod
    n::Int
    xn::Vector{T}
    yn::Matrix{T}
    c::Array{T,3}
    type::Symbol
end

"""
    InterpCubicSplines(x::AbstractVector, y::AbstractArray, type::Symbol=:Natural)

Construct a cubic spline interpolant from a set of data points `x` and their values `y`. 
Multi-dimensional splines can be constructed by passing `y` as a subtype of `AbstractMatrix`, 
such that each row contains a different set of values to be interpolated and the number of 
columns equals the number of data points.

Different boundary conditions can be applied based on the specified `type`: 

- **:Natural**: the second derivative of the first and the last polynomial are equal to 
    zero at the boundary points. 
- **:NotAKnot**: the third derivatives of the first and last two polynomials are 
    equal in the points where they meet each other.
- **:Periodic**: the first and second derivatives at the initial and final points are equal.
- **:Quadratic**: the first and the last polynomial are quadratic. 

"""
function InterpCubicSplines(x::AbstractVector, y::AbstractArray, type::Symbol=:Natural)

    # Check input validity 
    n = length(x)
    n < 4 && throw(ArgumentError("At least four points are needed."))

    ax = axes(y)
    N = length(ax)

    if N == 1
        ny = length(y)
    elseif N == 2
        ny = size(y, 2)
    else
        throw(ArgumentError("Arrays with more than 2 dimensions are not supported."))
    end

    ny != n && throw(ArgumentError("`x` and `y` must have the same length."))

    # Sort the arrays to guarantee that x is in ascending order
    idx = sortperm(x)

    xs = collect(x[idx])
    ys = reshape(collect(N == 1 ? y[idx] : y[:, idx]), N, n)

    # Compute the spline coefficients
    @views coeffs = vcat((_assemble_cspline(n, xs, ys[j, :], Val{type}())[:] for j in 1:N)...)
    T = eltype(coeffs)

    return InterpCubicSplines{T,N}(n, xs, ys, reshape(coeffs, (3, n - 1, N)), type)
end

Base.eltype(::InterpCubicSplines{T}) where {T} = T

""" 
    interpolate(sp::InterpCubicSplines, x::Number, flat::Bool=true)

Interpolate the cubic spline `sp` at point `x`. If the spline has a single 
dimension (e.g., `InterpCubicSpline{T, 1}`), a scalar value is returned. Otherwise,
an `SVector` is computed. 

If `x` is outside the boundary range of `sp` a flat extrapolation is used by default. 
If the `flat` argument is `false`, the first and last polynomials will be used to 
compute all the outside values.

"""
function jMath.interpolate(
    cs::InterpCubicSplines{T,1}, x::Number, flat::Bool=true
) where {T}
    @inbounds begin

        # Flat extrapolation settings
        if flat
            if x < cs.xn[1]
                return cs.yn[1, 1]
            elseif x > cs.xn[end]
                return cs.yn[1, end]
            end
        end

        # Search segment index
        j = max(2, min(searchsortedfirst(cs.xn, x), cs.n))

        # Horner polynomial evaluation
        δx = x - cs.xn[j - 1]
        @evalpoly δx cs.yn[1, j - 1] cs.c[1, j - 1] cs.c[2, j - 1] cs.c[3, j - 1]
    end
end

function jMath.interpolate(
    cs::InterpCubicSplines{T,N}, x::Number, flat::Bool=true
) where {T,N}
    @inbounds @views begin

        # Flat extrapolation settings
        @views if flat
            if x < cs.xn[1]
                return SVector{N}(cs.yn[:, 1])
            elseif x > cs.xn[end]
                return SVector{N}(cs.yn[:, end])
            end
        end

        # Search segment index
        j = max(2, min(searchsortedfirst(cs.xn, x), cs.n))
        δx = x - cs.xn[j - 1]

        # Horner polynomial evaluation
        return SVector{N}(
            cs.yn[i, j - 1] +
            δx * (cs.c[1, j - 1, i] + δx * (cs.c[2, j - 1, i] + δx * cs.c[3, j - 1, i])) for
            i in 1:N
        )
    end
end

# Compute cubic spline polynomial coefficients
function _assemble_cspline(n::Int, x::AbstractVector, y::AbstractVector, type)
    @inbounds @views begin
        δx = diff(x)
        δf = diff(y) ./ δx

        # Coefficient vector
        q = similar(δf, n)

        # Arrays storing diagonals, subdiagonal and superdiagonal
        dd = similar(δf, n)
        dl = similar(δf, n - 1)
        du = similar(δf, n - 1)

        # Compute the diagonal
        dd[2:(end - 1)] .= 2 * (δx[2:end] + δx[1:(end - 1)])

        # Compute the super and subdiagonals
        du[1] = 0
        du[2:end] .= δx[2:end]

        dl[1:(end - 1)] .= δx[1:(end - 1)]
        dl[end] = 0

        # Assemble bi coefficients vector
        q[2:(end - 1)] .= 3 * (δf[2:end] .- δf[1:(end - 1)])

        A = Tridiagonal(dl, dd, du)

        # Solve the system for the bi coefficients
        b = get_cspline_coefficients(type, A, q, δf, δx)

        # Polynomials coefficients matrix
        coeffs =
            hcat(
                δf .- δx .* (b[2:end] + 2 * b[1:(end - 1)]) / 3,
                b[1:(end - 1)],
                (b[2:end] - b[1:(end - 1)]) ./ (3 * δx),
            )'
    end

    return coeffs
end

# Apply natural boundary conditions, such that the second derivatives
# of the first and last polynomials are null at the start and end points.
@inbounds function get_cspline_coefficients(::Val{:Natural}, A, q, args...)
    A[1, 1] = 1
    A[end, end] = 1

    q[1] = 0.0
    q[end] = 0.0

    return A \ q
end

# Apply Not-A-Know boundary conditions, such that the third derivatives 
# where the first and last two polynomials match are equal
@inbounds function get_cspline_coefficients(::Val{:NotAKnot}, A, q, δf, δx)
    As = sparse(A)
    As[1, 1] = -δx[2]
    As[1, 2] = δx[1] + δx[2]
    As[1, 3] = -δx[1]

    As[end, end] = -δx[end - 1]
    As[end, end - 1] = δx[end] + δx[end - 1]
    As[end, end - 2] = -δx[end]

    q[1] = 0.0
    q[end] = 0.0

    return As \ q
end

# Apply periodic boundary conditions, such that the first and second derivatives
# at the first and last point are equal.
@inbounds function get_cspline_coefficients(::Val{:Periodic}, A, q, δf, δx, args...)
    As = sparse(A)
    As[1, 1] = 1
    As[1, end] = -1

    As[end, end] = 2δx[end]
    As[end, end - 1] = δx[end]
    As[end, 1] = 2δx[1]
    As[end, 2] = δx[1]

    q[1] = 0.0
    q[end] = 3 * (δf[1] - δf[end])

    return As \ q
end

# Apply quadratic boundary conditions, such that the first and last 
# polynomials are quadratic.
@inbounds function get_cspline_coefficients(::Val{:Quadratic}, A, q, args...)
    A[1, 1] = 1
    A[1, 2] = -1

    A[end, end] = 1
    A[end, end - 1] = -1

    q[1] = 0.0
    q[end] = 0.0

    return A \ q
end
