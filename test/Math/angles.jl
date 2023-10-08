
@testset "Angles" verbose = true begin
    # Test arcseconds to radians 
    @test arcsec2rad(3600 * 180) ≈ π atol = 1e-11 rtol = 1e-11
    @test rad2arcsec(arcsec2rad(3600 * 180)) ≈ 3600 * 180 atol = 1e-11 rtol = 1e-11
end;
