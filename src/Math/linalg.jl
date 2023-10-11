export skew, 
       unitvec, 
       cross3,
       unitcross, 
       projvec, projplane, 
       anglevec, anglevecd,
       angleplane, angleplaned

"""
    skew(a)
   
Create a skew matrix from the vector `v`.
"""
function skew(v::AbstractArray)
    @assert length(v) == 3
    @inbounds SMatrix{3, 3}(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0)
end


# -------------------------
# Vector normalisation

"""
    unitvec(v::AbstractVector)

Normalise the vector `v`.
"""
function unitvec(v::AbstractVector)
    @inbounds begin
        @fastmath r = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
        SA[v[1] / r, v[2] / r, v[3] / r]
    end
end

"""
    δunitvec(v::AbstractVector)

Compute the time derivative of a unit vector `v`.
"""
function δunitvec(v::AbstractVector)
    @inbounds begin
        r2 = v[1]^2 + v[2]^2 + v[3]^2
        @fastmath r = sqrt(r2)
        r3 = r2 * r

        δ = -(v[1] * v[4] + v[2] * v[5] + v[3] * v[6]) / r3

        SA[v[4] / r + δ * v[1], v[5] / r + δ * v[2], v[6] / r + δ * v[3]]
    end
end

"""
    δ²unitvec(v::AbstractVector)

Compute the 2nd-order time derivative of a unit vector `v`.
"""
function δ²unitvec(v::AbstractVector)
    @inbounds begin
        δ = v[1] * v[4] + v[2] * v[5] + v[3] * v[6]
        Δ = v[1] * v[7] + v[2] * v[8] + v[3] * v[9]

        r2 = v[1]^2 + v[2]^2 + v[3]^2
        @fastmath r = sqrt(r2)
        r3 = r2 * r

        dr² = v[4]^2 + v[5]^2 + v[6]^2

        A = 2 * δ / r3
        B = (dr² + Δ - 3δ^2 / r2) / r3

        a = v[7] / r - v[4] * A - v[1] * B
        b = v[8] / r - v[5] * A - v[2] * B
        c = v[9] / r - v[6] * A - v[3] * B
    end

    return SA[a, b, c]
end

"""
    δ³unitvec(v::AbstractVector)

Compute the 3rd-order time derivative of a unit vector `v`.
"""
function δ³unitvec(v::AbstractVector)
    @inbounds begin
        δ = v[1] * v[4] + v[2] * v[5] + v[3] * v[6]
        Δ = v[1] * v[7] + v[2] * v[8] + v[3] * v[9]
        θ = v[4] * v[7] + v[5] * v[8] + v[6] * v[9]
        φ = v[1] * v[10] + v[2] * v[11] + v[3] * v[12]

        r2 = v[1]^2 + v[2]^2 + v[3]^2
        @fastmath r = sqrt(r2)
        r3 = r2 * r

        dr² = v[4]^2 + v[5]^2 + v[6]^2
        δ² = δ^2

        A = 3δ / r3
        B = 3 * (Δ + dr² - 3δ² / r2) / r3
        C = (3θ + φ - 3δ / r2 * (3Δ + 3dr² - 5δ² / r2)) / r3

        a = v[10] / r - v[7] * A - v[4] * B - v[1] * C
        b = v[11] / r - v[8] * A - v[5] * B - v[2] * C
        c = v[12] / r - v[9] * A - v[6] * B - v[3] * C
    end

    return SA[a, b, c]
end


# -------------------------
# Cross products 

"""

    cross3(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y`, considering only their first 3 elements.
"""
function cross3(x::AbstractVector, y::AbstractVector)
    @inbounds begin
        u = x[2] * y[3] - x[3] * y[2]
        v = x[3] * y[1] - y[3] * x[1]
        w = x[1] * y[2] - x[2] * y[1]
    end

    return SA[u, v, w]
end

"""
    cross6(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y` and its time derivative. 

!!! note
    `x` and `y` must be 6-elements state vectors, with the last elements of each vector 
    representing the time derivatives of the first three.
"""
function cross6(x::AbstractVector, y::AbstractVector)
    @inbounds begin
        u = x[2] * y[3] - x[3] * y[2]
        v = x[3] * y[1] - y[3] * x[1]
        w = x[1] * y[2] - x[2] * y[1]

        δu = x[5] * y[3] + x[2] * y[6] - x[6] * y[2] - x[3] * y[5]
        δv = x[6] * y[1] + x[3] * y[4] - x[4] * y[3] - x[1] * y[6]
        δw = x[4] * y[2] + x[1] * y[5] - x[5] * y[1] - x[2] * y[4]
    end

    return SA[u, v, w, δu, δv, δw]
end

"""
    cross9(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y` and its 1st and 2nd-order time derivatives. 

!!! note
    `x` and `y` must be 9-elements state vectors (position, velocity and acceleration)
"""
function cross9(x::AbstractVector, y::AbstractVector)
    @inbounds begin
        u = x[2] * y[3] - x[3] * y[2]
        v = x[3] * y[1] - y[3] * x[1]
        w = x[1] * y[2] - x[2] * y[1]

        δu = x[5] * y[3] + x[2] * y[6] - x[6] * y[2] - x[3] * y[5]
        δv = x[6] * y[1] + x[3] * y[4] - x[4] * y[3] - x[1] * y[6]
        δw = x[4] * y[2] + x[1] * y[5] - x[5] * y[1] - x[2] * y[4]

        δ²u =
            x[8] * y[3] + 2x[5] * y[6] + x[2] * y[9] - x[9] * y[2] - 2x[6] * y[5] -
            x[3] * y[8]
        δ²v =
            x[9] * y[1] + 2x[6] * y[4] + x[3] * y[7] - x[7] * y[3] - 2x[4] * y[6] -
            x[1] * y[9]
        δ²w =
            x[7] * y[2] + 2x[4] * y[5] + x[1] * y[8] - x[8] * y[1] - 2x[5] * y[4] -
            x[2] * y[7]
    end

    return SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w]
end

"""
    cross12(x::AbstractVector, y::AbstractVector)

Compute the cross product between `x` and `y` and its 1st, 2nd and 3rd order time derivatives. 

!!! note
    `x` and `y` must be 12-elements state vectors (position, velocity and acceleration)
"""
function cross12(x::AbstractVector, y::AbstractVector)
    @inbounds begin
        u = x[2] * y[3] - x[3] * y[2]
        v = x[3] * y[1] - y[3] * x[1]
        w = x[1] * y[2] - x[2] * y[1]

        δu = x[5] * y[3] + x[2] * y[6] - x[6] * y[2] - x[3] * y[5]
        δv = x[6] * y[1] + x[3] * y[4] - x[4] * y[3] - x[1] * y[6]
        δw = x[4] * y[2] + x[1] * y[5] - x[5] * y[1] - x[2] * y[4]

        δ²u =
            x[8] * y[3] + 2x[5] * y[6] + x[2] * y[9] - x[9] * y[2] - 2x[6] * y[5] -
            x[3] * y[8]
        δ²v =
            x[9] * y[1] + 2x[6] * y[4] + x[3] * y[7] - x[7] * y[3] - 2x[4] * y[6] -
            x[1] * y[9]
        δ²w =
            x[7] * y[2] + 2x[4] * y[5] + x[1] * y[8] - x[8] * y[1] - 2x[5] * y[4] -
            x[2] * y[7]

        δ³u =
            x[11] * y[3] + 3x[8] * y[6] + 3x[5] * y[9] + x[2] * y[12] + -x[12] * y[2] -
            3x[9] * y[5] - 3x[6] * y[8] - x[3] * y[11]

        δ³v =
            x[12] * y[1] + 3x[9] * y[4] + 3x[6] * y[7] + x[3] * y[10] + -x[10] * y[3] -
            3x[7] * y[6] - 3x[4] * y[9] - x[1] * y[12]

        δ³w =
            x[10] * y[2] + 3x[7] * y[5] + 3x[4] * y[8] + x[1] * y[11] + -x[11] * y[1] -
            3x[8] * y[4] - 3x[5] * y[7] - x[2] * y[10]
    end

    return SA[u, v, w, δu, δv, δw, δ²u, δ²v, δ²w, δ³u, δ³v, δ³w]
end


# -------------------------
# Other algebra routines


"""
    unitcross(v1::AbstractArray, v2::AbstractArray)

Compute the normalized cross product of `v1` and `v2`.
"""
function unitcross(v1::AbstractArray, v2::AbstractArray)
    tmp = cross(v1, v2)
    return unitvec(tmp)
end

"""
    projvec(v1::AbstractArray, v2::AbstractArray)

Compute the orthogonal projection of a vector `v1` on a vector `v2`.
"""
function projvec(v1::AbstractArray, v2::AbstractArray)
    v2n2 = dot(v2, v2)
    if v2n2 == 0.
        return SVector{3}(0., 0., 0.)
    end
    return dot(v1, v2) * v2/v2n2
end

"""
    projplane(v1::AbstractArray, n::AbstractArray)

Compute the orthogonal projection of a vector `v1` on a plane with normal `n`.
"""
function projplane(v1::AbstractArray, n::AbstractArray)
    return v1 - projvec(v1, n)
end

"""
    anglevec(v1::AbstractArray, v2::AbstractArray)

Compute the angle between two vectors `v1` and `v2`, in rad.
"""
function anglevec(v1::AbstractArray, v2::AbstractArray)
    tmp = dot(unitvec(v1), unitvec(v2))
    return mod(acos(tmp), π)
end

"""
    anglevecd(v1::AbstractArray, v2::AbstractArray)

Compute the angle between two vectors `v1` and `v2`, in deg.
"""
function anglevecd(v1::AbstractArray, v2::AbstractArray)
    return rad2deg(anglevec(v1, v2))
end

"""
    angleplane(v1::AbstractArray, n::AbstractArray)

Compute the angle between a vector `v1` and (its projection on) a plane with
normal `n`, in rad. 
"""
function angleplane(v1::AbstractArray, n::AbstractArray)
    return anglevec(projplane(v1, n), v1)
end

"""
    angleplaned(v1::AbstractArray, n::AbstractArray)

Compute the angle between a vector `v1` and (its projection on) a plane with
normal `n`, in deg. 
"""
function angleplaned(v1::AbstractArray, n::AbstractArray)
    return rad2deg(angleplane(v1, n))
end