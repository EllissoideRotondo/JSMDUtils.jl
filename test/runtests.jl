using JSMDUtils
using Test 

@testset "JSMDUtils" verbose = true begin
    @eval begin
        modules = [:Math]
        for m in modules
            @testset "$m" verbose = true begin
                include("$m/$m.jl")
            end
        end
    end
end