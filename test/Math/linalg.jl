
# Define two non parallel vectors and their 1st, 2nd and 3rd order time derivatives! 
function get_vector(t)
    return SA[
        cos(3t),
        t * sin(t),
        t^2 * cos(t),
        -3sin(3t),
        sin(t) + t * cos(t),
        2t * cos(t) - t^2 * sin(t),
        -9cos(3t),
        2cos(t) - t * sin(t),
        2cos(t) - 4t * sin(t) - t^2 * cos(t),
        27sin(3t),
        -3sin(t) - t * cos(t),
        -6sin(t) - 6t * cos(t) + t^2 * sin(t),
    ]
end

get_vector2(t) = SA[t^3, t^2, t, 3t^2, 2t, 1, 6t, 2, 0, 6, 0, 0]

# Function to compute derivatives of normalization 
δnorm(t) = get_vector(t)[1:3] / norm(get_vector(t)[1:3])

# Function to compute derivatives of cross product 
δcross(t) = cross(get_vector(t)[1:3], get_vector2(t)[1:3])

@testset "Skew-matrix" begin
    
    S = Math.skew([1,2,3])
    @test sum(diag(S)) == 0
    @test S[1, 2] == -3
    @test S[1, 3] == 2 
    @test S[2, 1] == 3
    @test S[2, 3] == -1
    @test S[3, 1] == -2 
    @test S[3, 2] == 1

end


# Normalisation\cross product routines 
atol = 1e-11

# Unit vector derivatives 
# -----------------------
@testset "Normalisation" verbose = true begin

    # Check vector normalisation
    a, b = [-2.0, 3.0, 1.0, 4], SA[1, 2, 3]
    @test Math.unitvec(b) ≈ b / norm(b) atol = atol
    @test Math.unitvec(a) ≈ SA[a[1:3]...] / norm(a[1:3]) atol = atol

    for _ in 1:100
        θ = rand()
        # Compute 1st, 2nd and 3rd order derivatives of unit vector 
    
        v = get_vector(θ)
        @test Math.δunitvec(v) ≈ D¹(δnorm, θ) atol = atol
        @test Math.δ²unitvec(v) ≈ D²(δnorm, θ) atol = atol
        @test Math.δ³unitvec(v) ≈ D³(δnorm, θ) atol = atol
    end

end


# Cross product derivatives 
# -------------------------  
@testset "Cross Product" verbose = true begin

    for _ in 1:100
        
        θ = rand()
        a, b = get_vector(θ), get_vector2(θ)

        for (i, fcn) in enumerate([Math.cross6, Math.cross9, Math.cross12])
            cp = fcn(a, b)

            @test cp[1:3] ≈ cross(a[1:3], b[1:3]) atol = atol # Cross product 
            @test cp[4:6] ≈ D¹(δcross, θ) atol = atol # 1st order derivative   

            if i > 1 # 2nd order derivative 
                @test cp[7:9] ≈ D²(δcross, θ) atol = atol
            end

            if i > 2 # 3rd order derivative 
                @test cp[10:12] ≈ D³(δcross, θ) atol = atol
            end
        end 
    end
end


@testset "projvec" begin 
    I = [1., 0., 0.]
    J = [0., 1., 0.]
    @test Math.projvec(I, J) == zeros(3)

    v = [cosd(30), sind(30), 0.]
    @test Math.projvec(v, I) == [v[1], 0., 0.]
    @test Math.projvec(v, J) == [0., v[2], 0.]

    n = [0., 0., 1.]
    @test Math.projplane(I, n) == I
    @test Math.projplane(J, n) == J
    @test Math.projplane(v, n) == v

    v = [cosd(30)*cosd(45), sind(30)*cosd(45), sind(45)]
    @test Math.projplane(v, n) == [v[1], v[2], 0.]
end

@testset "anglevec/anglevecd" begin
    I = [1., 0., 0.]
    J = [0., 1., 0.]
    @test Math.anglevec(I, J) == π/2
    @test Math.anglevec(I, -J) == π/2
    @test Math.anglevec(-I, -J) == π/2

    v = [cosd(30), sind(30), 0.]
    @test Math.anglevec(I, v) ≈ π/6
    @test Math.anglevec(v, I) ≈ π/6
    @test Math.anglevecd(J, v) ≈ 60.
end

@testset "angleplane/angleplaned" begin
    
    v = [cosd(30)*cosd(45), sind(30)*cosd(45), sind(45)]
    @test Math.angleplane(v, [0., 0., 1.]) ≈ π/4
    @test Math.angleplaned(v, [0., 0., 1.]) ≈ 45.
    
end

