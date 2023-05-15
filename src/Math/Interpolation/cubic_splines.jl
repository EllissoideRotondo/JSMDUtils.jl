export InterpCubicSplines

struct InterpCubicSplines{T,N} <: jMath.AbstractInterpolationMethod
    n::Int
    xn::Vector{T}
    yn::Matrix{T}
    c::Array{T,3}
    type::Symbol
end

"""
    InterpCubicSplines(x::AbstractVector, y::AbstractArray, type::Symbol)
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
    @views coeffs = vcat((_assemble_cspline(n, xs, ys[j, :], Val(type))[:] for j in 1:N)...)
    T = eltype(coeffs)

    return InterpCubicSplines{T,N}(n, xs, ys, reshape(coeffs, (4, n - 1, N)), type)
end

""" 
    interpolate(::InterpCubicSplines, x::Number)
"""
function jMath.interpolate(
    cs::InterpCubicSplines{T,1}, x::Number, flat::Bool=true
) where {T}

    # Flat extrapolation settings
    if flat
        if x < cs.xn[1]
            return cs.yn[1, 1]
        elseif x > cs.xn[end]
            return cs.yn[1, end]
        end
    end

    @inbounds begin
        j = max(2, min(searchsortedfirst(cs.xn, x), cs.n))

        # Horner polynomial evaluation
        δx = x - cs.xn[j - 1]
        return cs.c[1, j - 1] +
               δx * (cs.c[2, j - 1] + δx * (cs.c[3, j - 1] + δx * cs.c[4, j - 1]))
    end
end

function jMath.interpolate(
    cs::InterpCubicSplines{T,N}, x::Number, flat::Bool=true
) where {T,N}

    # Flat extrapolation settings
    if flat
        if x < cs.xn[1]
            return SVector{N}(yi for yi in @view cs.yn[1:end, 1])
        elseif x > cs.xn[end]
            return SVector{N}(yi for yi in @view cs.yn[1:end, end])
        end
    end

    @inbounds begin
        j = max(2, min(searchsortedfirst(cs.xn, x), cs.n))

        # Horner polynomial evaluation
        δx = x - cs.xn[j - 1]
        return SVector{N}(
            cs.c[1, j - 1, i] +
            δx * (cs.c[2, j - 1, i] + δx * (cs.c[3, j - 1, i] + δx * cs.c[4, j - 1, i])) for
            i in 1:N
        )
    end
end

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
                y[1:(end - 1)],
                δf .- δx .* (b[2:end] + 2 * b[1:(end - 1)]) / 3,
                b[1:(end - 1)],
                (b[2:end] - b[1:(end - 1)]) ./ (3 * δx),
            )'
    end

    return coeffs
end

# Apply natural boundary conditions, such that the second derivatives
# of the first and last polynomials are null at the start and end points.
function get_cspline_coefficients(::Val{:Natural}, A, q, args...)
    A[1, 1] = 1
    A[end, end] = 1

    q[1] = 0.0
    q[end] = 0.0

    return A \ q
end

# Apply Not-A-Know boundary conditions, such that the third derivatives 
# where the first and last two polynomials match are equal
function get_cspline_coefficients(::Val{:NotAKnot}, A, q, δf, δx)
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
function get_cspline_coefficients(::Val{:Periodic}, A, q, δf, δx, args...)
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
function get_cspline_coefficients(::Val{:Quadratic}, A, q, args...)
    A[1, 1] = 1
    A[1, 2] = -1

    A[end, end] = 1
    A[end, end - 1] = -1

    q[1] = 0.0
    q[end] = 0.0

    return A \ q
end
