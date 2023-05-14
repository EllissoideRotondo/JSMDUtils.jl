using JSMDUtils
using Test

import JSMDInterfaces.Ephemeris as jEphem

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
end;
