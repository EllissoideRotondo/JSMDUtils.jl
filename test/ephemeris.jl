
@testset "Ephemeris" verbose = true begin
    np = JSMDUtils.NullEphemerisProvider()

    @test jEphem.ephem_timescale(np) == 1
end
