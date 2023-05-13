
@testset "Ephemeris" verbose=true begin 

    np = JSMDUtils.NullEphemerisProvider()

    # TODO: with JSMDInterfaces change to that caller
    @test JSMDUtils.ephem_timescale(np) == 1

end