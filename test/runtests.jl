using JSMDUtils
using Test

using ForwardDiff
using LazyArtifacts
using LinearAlgebra
using ReferenceFrameRotations
using StaticArrays

using JSMDUtils.Autodiff
using JSMDUtils.Math: Math, D¹, D², D³, arcsec2rad, rad2arcsec, 
                      angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm

import JSMDInterfaces.Ephemeris as jEphem
import JSMDInterfaces.FilesIO as jIO
import JSMDInterfaces.Math as jMath

import JSMDUtils.FileUtils

@testset "Download all artifacts" begin
    @info artifact"testdata"
    @info "All artifacts downloaded"
end

@testset "JSMDUtils" verbose = true begin
    @eval begin
        modules = [:Math]
        for m in modules
            @testset "$m" verbose = true begin
                include("$m/$m.jl")
            end
        end
    end

    include("ephemeris.jl")
    include("format.jl")
    include("FileUtils.jl")
    include("Autodiff.jl")
    
end;
