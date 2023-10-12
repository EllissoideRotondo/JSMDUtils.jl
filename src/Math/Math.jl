module Math

using LinearAlgebra
using ForwardDiff: derivative, ForwardDiff
using PrecompileTools
using ReferenceFrameRotations: DCM
using StaticArrays
using SparseArrays

import JSMDInterfaces.Math as jMath

include("angles.jl")
include("derivatives.jl")
include("linalg.jl")
include("rotations.jl")

include(joinpath("Interpolation", "Interpolation.jl"))

# Precompilation routines 
PrecompileTools.@setup_workload begin

    x3 = rand(3)
    x6 = rand(6)
    x9 = rand(9)
    x12 = rand(12)

    x3s = SA[rand(3)...]
    x6s = SA[rand(6)...]
    x9s = SA[rand(9)...]
    x12s = SA[rand(12)...]

    vfcns = (
        (unitvec, cross3),
        (δunitvec, cross6),
        (δ²unitvec, cross9),
        (δ³unitvec, cross12),
    )

    vects = ((x3, x3s), (x6, x6s), (x9, x9s), (x12, x12s))

    PrecompileTools.@compile_workload begin

        # Precompile angles routines 
        arcsec2rad(0)
        rad2arcsec(0)
        
        arcsec2rad(0.0)
        rad2arcsec(0.0)

        # Precompile vector routines
        for (vfcn, vx) in zip(vfcns, vects)
            for x in vx
                vfcn[1](x)
                vfcn[2](x, x)
            end
        end

        # Precompile rotation routines 
        skew(x3)
        skew(x3s)

        angle_to_δdcm(x3, :Z)
        angle_to_δdcm(x3, x3, :ZX)
        angle_to_δdcm(x3, x3, x3, :ZXZ)

        angle_to_δ²dcm(x3, :Z)
        angle_to_δ²dcm(x3, x3, :ZX)
        angle_to_δ²dcm(x3, x3, x3, :ZXZ)

        angle_to_δ³dcm(x6, :Z)
        angle_to_δ³dcm(x6, x6, :ZX)
        angle_to_δ³dcm(x6, x6, x6, :ZXZ)

        _3angles_to_δdcm(x6, :ZXZ)
        _3angles_to_δdcm(x6s, :ZXZ)

        _3angles_to_δ²dcm(x9, :ZXZ)
        _3angles_to_δ²dcm(x9s, :ZXZ)

        _3angles_to_δ³dcm(x12, :ZXZ)
        _3angles_to_δ³dcm(x12s, :ZXZ)

    end

end


end
