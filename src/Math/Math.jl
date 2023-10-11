module Math

using LinearAlgebra
using ForwardDiff: derivative, ForwardDiff
using ReferenceFrameRotations: DCM
using StaticArrays
using SparseArrays

import JSMDInterfaces.Math as jMath

include("angles.jl")
include("derivatives.jl")
include("linalg.jl")
include("rotations.jl")

include(joinpath("Interpolation", "Interpolation.jl"))

end
