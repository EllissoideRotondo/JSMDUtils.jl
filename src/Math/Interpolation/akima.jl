export InterpAkima

"""
    InterpAkima{T, N} <: AbstractInterpolationMethod

Type storing an Akima spline nodes and coefficient. `T` is the interpolation data type 
and `N` is the spline dimension (i.e., the number of interpolated functions).

### Fields 
- `n` -- Number of node points.
- `xn` -- Interpolated node points. 
- `yn` -- Node points function values 
- `c` -- Akima polynomial coefficients. 
"""
struct InterpAkima{T,N} <: jMath.AbstractInterpolationMethod
    n::Int
    xn::Vector{T}
    yn::Matrix{T}
    c::Array{T,3}
end

"""
    InterpAkima(x::AbstractVector, y::AbstractArray)

Construct an Akima spline interpolant from a set of data points `x` and their values `y`. 
Multi-dimensional splines can be constructed passing `y` as a subtype of `AbstractMatrix` such 
that each row contains a different set of values to be interpolated and the number of 
columns equals the number of data points.

### References
-  Akima, H. (1970),  A New Method of Interpolation and Smooth Curve Fitting Based on Local
    Procedures, Journal of the ACM, [DOI:](https://dl.acm.org/doi/10.1145/321607.321609) 
"""
function InterpAkima(x::AbstractVector, y::AbstractArray)

    # Check input validity 
    n = length(x)

    if n < 5
        throw(ArgumentError("At least 5 points are needed to construct an Akima spline."))
    end

    ax = axes(y)
    N = length(ax)

    if N == 1
        ny = length(y)
    elseif N == 2
        ny = size(y, 2)
    else
        throw(ArgumentError("Arrays with more than 2 dimensions are not supported"))
    end

    ny != n && throw(ArgumentError("`x` and `y` must have the same length."))

    # Sort the arrays to guarantee that x is in ascending order 
    idx = sortperm(x)

    xs = collect(x[idx])
    ys = reshape(collect(N == 1 ? y[idx] : y[:, idx]), N, n)

    # Compute the spline coefficients 
    @views coeffs = vcat((_assemble_akima(n, xs, ys[j, :])[:] for j in 1:N)...)
    T = eltype(coeffs)

    return InterpAkima{T,N}(n, xs, ys, reshape(coeffs, (3, n - 1, N)))
end

Base.eltype(::InterpAkima{T}) where {T} = T

"""
    interpolate(ak::InterpAkima, x::Number, flat::Bool=true)

Interpolate the Akima spline `ak` at point `x`. If the spline has a single 
dimension (e.g., `InterpAkima{T, 1}`), a scalar value is returned. Otherwise,
an `SVector` is computed. 

If `x` is outside the boundary range of `sp` a flat extrapolation is used by default. 
If the `flat` argument is `false`, the first and last polynomials will be used to 
compute all the outside values.
"""
function jMath.interpolate(ak::InterpAkima{T,1}, x::Number, flat::Bool=true) where {T}
    @inbounds begin

        # Flat extrapolation settings 
        if flat
            if x < ak.xn[1]
                return ak.yn[1, 1]
            elseif x > ak.xn[end]
                return ak.yn[1, end]
            end
        end

        # Search segment index
        j = max(2, min(searchsortedfirst(ak.xn, x), ak.n))

        # Horner polynomial evaluation 
        δx = x - ak.xn[j - 1]
        @evalpoly δx ak.yn[1, j - 1] ak.c[1, j - 1] ak.c[2, j - 1] ak.c[3, j - 1]
    end
end

function jMath.interpolate(ak::InterpAkima{T,N}, x::Number, flat::Bool=true) where {T,N}
    @inbounds @views begin

        # Flat extrapolation settings 
        if flat
            if x < ak.xn[1]
                return SVector{N}(ak.yn[:, 1])
            elseif x > ak.xn[end]
                return SVector{N}(ak.yn[:, end])
            end
        end

        # Search segment index
        j = max(2, min(searchsortedfirst(ak.xn, x), ak.n))
        δx = x - ak.xn[j - 1]

        # Horner polynomial evaluation 
        return SVector{N}(
            ak.yn[i, j - 1] +
            δx * (ak.c[1, j - 1, i] + δx * (ak.c[2, j - 1, i] + δx * ak.c[3, j - 1, i])) for
            i in 1:N
        )
    end
end

# Compute akima polynomial coefficients
@inbounds @views function _assemble_akima(n::Int, x::AbstractVector, y::AbstractVector)
    δx = diff(x)
    δf = diff(y) ./ δx

    # Compute slopes of each line segment 
    m = similar(δf, n + 3)
    m[3:(end - 2)] = δf

    # Compute the slopes of the first and last line segments (with quadratic polynomial)
    m[2] = 2m[3] - m[4]
    m[1] = 2m[2] - m[3]

    m[end - 1] = 2m[end - 2] - m[end - 3]
    m[end] = 2m[end - 1] - m[end - 2]

    # If m1 == m2 != m3 == m4, the slope at each node is not defined.
    # In such cases this value will be used 
    s = (m[3:(end - 1)] + m[2:(end - 2)]) / 2

    # Compute the slope denominator:
    δm = abs.(diff(m))

    f1 = δm[3:(n + 2)]
    f2 = δm[1:n]

    f12 = f1 + f2

    # Retrieve the node indexes where the slope is not defined: 
    idx = findall(f12 .> 0)

    # Compute the slope of these points
    s[idx] = (f1[idx] .* m[idx .+ 1] .+ f2[idx] .* m[idx .+ 2]) ./ f12[idx]

    # Compute the remaining coefficients 
    coeffs =
        hcat(
            s[1:(end - 1)],
            (3m[3:(end - 2)] .- 2s[1:(end - 1)] .- s[2:end]) ./ δx,
            (s[1:(end - 1)] .+ s[2:end] .- 2m[3:(end - 2)]) ./ δx .^ 2,
        )'

    return coeffs
end
