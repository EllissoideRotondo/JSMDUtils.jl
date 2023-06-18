using JSMDUtils.Math
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

@testset "unitvec" begin
    I = [1., 0., 0.]
    J = [0., 1., 0.]
    K = [0., 0., 1.]
    @test unitvec(I) == I
    @test unitvec(J) == J
    @test unitvec(K) == K

    v = rand(3)
    vn = norm(v)
    @test unitvec(v) == v/vn
end

@testset "unitcross" begin 
    I = [1., 0., 0.]
    J = [0., 1., 0.]
    K = [0., 0., 1.]
    @test unitcross(I, J) == K

    v1 = [2., 0., 0.]
    @test unitcross(v1, J) == K 
end

@testset "projvec" begin 
    I = [1., 0., 0.]
    J = [0., 1., 0.]
    @test projvec(I, J) == zeros(3)

    v = [cosd(30), sind(30), 0.]
    @test projvec(v, I) == [v[1], 0., 0.]
    @test projvec(v, J) == [0., v[2], 0.]

    n = [0., 0., 1.]
    @test projplane(I, n) == I
    @test projplane(J, n) == J
    @test projplane(v, n) == v

    v = [cosd(30)*cosd(45), sind(30)*cosd(45), sind(45)]
    @test projplane(v, n) == [v[1], v[2], 0.]
end

@testset "anglevec/anglevecd" begin
    I = [1., 0., 0.]
    J = [0., 1., 0.]
    @test anglevec(I, J) == π/2
    @test anglevec(I, -J) == π/2
    @test anglevec(-I, -J) == π/2

    v = [cosd(30), sind(30), 0.]
    @test anglevec(I, v) ≈ π/6
    @test anglevec(v, I) ≈ π/6
    @test anglevecd(J, v) ≈ 60.
end

@testset "angleplane/angleplaned" begin
    
    v = [cosd(30)*cosd(45), sind(30)*cosd(45), sind(45)]
    @test angleplane(v, [0., 0., 1.]) ≈ π/4
    @test angleplaned(v, [0., 0., 1.]) ≈ 45.
    
end