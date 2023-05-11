
# TODO: remove when JSMDInterfaces becomes available
abstract type AbstractEphemerisProvider end 

"""
    EmptyEphemerisProvider <: AbstractEphemerisProvider

Empty provider to initialise the frame system without loading 
ephemeris files. 
"""
struct NullEphemerisProvider <: AbstractEphemerisProvider end

ephem_timescale(::NullEphemerisProvider) = TDB
