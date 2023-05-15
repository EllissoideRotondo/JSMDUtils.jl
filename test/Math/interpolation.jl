
@testset "CubicSplines" verbose = true begin
    xn = LinRange(0, 2 * π, 10)
    yn = sin.(xn)

    for type in (:Natural, :Quadratic, :Periodic, :NotAKnot)
        cs = JSMDUtils.Math.InterpCubicSplines(xn, yn, type)
        @test cs.type == type

        # Test node values
        for j in eachindex(xn)
            @test jMath.interpolate(cs, xn[j]) ≈ yn[j] atol = 1e-11 rtol = 1e-11
        end

        df = x -> jMath.interpolate(cs, x)

        # Test 1st derivatives to be continuous
        for j in 2:(length(xn) - 1)
            @test D¹(df, xn[j] - 1e-15) ≈ D¹(df, xn[j] + 1e-15) atol = 1e-11 rtol = 1e-11
        end

        # Test 2nd derivatives to be continuous
        for j in 2:(length(xn) - 1)
            @test D²(df, xn[j] - 1e-15) ≈ D²(df, xn[j] + 1e-15) atol = 1e-11 rtol = 1e-11
        end

        if type == :Natural
            # 2nd derivative at first and last point must be null 
            @test D²(df, xn[1]) ≈ 0 atol = 1e-11 rtol = 1e-11

        elseif type == :Quadratic
            # First and last polynomials are quadratic, so third derivativse must be null
            @test D³(df, rand([xn[1], xn[2]])) ≈ 0 atol = 1e-11 rtol = 1e-11
            @test D³(df, rand([xn[end - 1], xn[end]])) ≈ 0 atol = 1e-11 rtol = 1e-11

        elseif type == :Periodic
            # 1st and 2nd derivatives at first and last point must be equal
            @test D¹(df, xn[1]) ≈ D¹(df, xn[end]) atol = 1e-11 rtol = 1e-11
            @test D²(df, xn[1]) ≈ D²(df, xn[end]) atol = 1e-11 rtol = 1e-11

        elseif type == :NotAKnot
            # Third derivatives between first two and last two must be equal 
            @test D³(df, xn[2] - 1e-15) ≈ D³(df, xn[2] + 1e-15) atol = 1e-11 rtol = 1e-11
            @test D³(df, xn[end - 1] - 1e-15) ≈ D³(df, xn[end - 1] + 1e-15) atol = 1e-11 rtol =
                1e-11
        end
    end
end
