module Math 

using ForwardDiff
using ForwardDiff: derivative

include("angles.jl")
include("derivatives.jl")

include(joinpath("Interpolation", "Interpolation.jl"))

end