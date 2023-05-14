
"""
    EmptyEphemerisProvider <: AbstractEphemerisProvider

Empty provider to initialise the frame system without loading 
ephemeris files. 
"""
struct NullEphemerisProvider <: jEphem.AbstractEphemerisProvider end

# Sets default timescale for NullEphemerisProvider to TDB
jEphem.ephem_timescale(::NullEphemerisProvider) = 1
