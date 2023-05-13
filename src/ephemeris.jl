
# TODO: remove when JSMDInterfaces becomes available
abstract type AbstractEphemerisProvider end 

"""
    EmptyEphemerisProvider <: AbstractEphemerisProvider

Empty provider to initialise the frame system without loading 
ephemeris files. 
"""
struct NullEphemerisProvider <: AbstractEphemerisProvider end

# Sets default timescale for NullEphemerisProvider to TDB
ephem_timescale(::NullEphemerisProvider) = 1
