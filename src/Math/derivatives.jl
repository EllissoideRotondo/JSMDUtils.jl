
"""
    D¹(f, x::Real)

Return `df/dx` evaluated at `x` using ForwardDiff, assuming `f` is called as `f(x)`.

This method assumes that `isa(f(x), Union{Real, AbstractArray})`.
"""
@inline D¹(f, t) = derivative(f, t)

"""
    D²(f, x::Real)

Return `d²f/dx²` evaluated at `x` using ForwardDiff, assuming `f` is called as `f(x)`.

This method assumes that `isa(f(x), Union{Real, AbstractArray})`.
"""
@inline D²(f, t) = derivative(τ -> derivative(f, τ), t)

"""
    D³(f, x::Real)

Return `d³f/dx³` evaluated at `x` using ForwardDiff, assuming `f` is called as `f(x)`.

This method assumes that `isa(f(x), Union{Real, AbstractArray})`.
"""
@inline D³(f, t) = derivative(κ -> derivative(τ -> derivative(f, τ), κ), t)

# Overload of extract_derivative to allow AD when output types are tuples 
# rather than subtypes of AbstractArrays
@inline function ForwardDiff.extract_derivative(::Type{T}, y::Tuple) where {T}
    return map(d -> ForwardDiff.extract_derivative(T, d), y)
end
