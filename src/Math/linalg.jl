export skew

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