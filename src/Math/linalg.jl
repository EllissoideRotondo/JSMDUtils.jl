export skew, 
       unitvec, 
       unitcross, 
       projvec, projplane, 
       anglevec, anglevecd,
       angleplane, angleplaned

@inbounds function skew(v::AbstractArray{N}) where {N<:Number}
    @assert length(v) == 3
    return SMatrix{3, 3, N, 9}(
          0.,    v[3],   -v[2], 
       -v[3],      0.,    v[1],
        v[2],   -v[1],      0.
    )
end

"""
    skew(v::AbstractArray{N}) where {N<:Number}

Compute a skew-symmetric matrix from the vector `v`. 
"""
skew

"""
    unitvec(v::AbstractArray)

Compute the normalized vector of `v`.
"""
function unitvec(v::AbstractArray)
    nv = norm(v)
    if nv == 1.
        return v
    else 
        return v/nv
    end
end

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