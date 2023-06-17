using JSMDUtils.Math: skew
using LinearAlgebra

@testset "Skew-matrix" begin
    
    S = skew([1,2,3])
    @test sum(diag(S)) == 0
    @test S[1, 2] == -3
    @test S[1, 3] == 2 
    @test S[2, 1] == 3
    @test S[2, 3] == -1
    @test S[3, 1] == -2 
    @test S[3, 2] == 1

end