
abstract type AbstractProbe end

abstract type AbstractProbe1d <: AbstractProbe end
abstract type AbstractProbe2d <: AbstractProbe end

abstract type AbstractProbe3d <: AbstractProbe end


mutable struct Probe1d{Anem<:AbstractThermalAnemometer,Setup} <: AbstractProbe1d
    "Anemometer sensor information"
    sensor::Anem
    "Anemometer setup"
    setup::Setup
end

Base.broadcastable(sensor::Probe1d) = Ref(sensor)

"""
    `Probe1d(sensor, cal, setup)`

Creates a 1d probe. This completely characterizes a 1d probe. It includes
information on calibration, bridge configuration, output filter, etc.

# Arguments:

 - `sensor`: an `AbstractThermalAnemometer` object
 - `cal`: `CalibrCurve` object

"""
#function Probe1d(sensor, cal, setup="")
#    return Probe1d(sensor, cal, setup)
#end


resistor(w::Probe1d) = resistor(w.sensor)
sensor(w::Probe1d) = w.sensor
gain(w::Probe1d) = gain(w.sensor)
reftemp(w::Probe1d) = reftemp(w.sensor)
resistance(w::Probe1d) = resistance(w.sensor)

# Just get stuff for CTA
temperature(w::Probe1d{<:AbstractCTA}) = temperature(w.sensor)
refresist(w::Probe1d{<:AbstractCTA}) = refresist(w.sensor)
overheatratio(w::Probe1d{<:AbstractCTA}) = overheatratio(w.sensor)
overheatratio(w::Probe1d{<:AbstractCTA},T) =
    overheatratio(w.sensor, T)

overtemp(w::Probe1d{<:AbstractCTA}) = overtemp(w.sensor)
overtemp(w::Probe1d{<:AbstractCTA}, T) = overtemp(w.sensor, T)

# Stuff for CCA
current(w::Probe1d{<:AbstractCCA}) = current(w.sensor)


# Now let's actually do something with `Probe1d`

correct(w::Probe1d, E::Real; kw...) = correct(sensor(w), args...; kw...)
velocity(w::Probe1d, E::Real; kw...) = velocity(sensor(w), E; kw...)
velocity(w::Probe1d, E::Real, fc::CorrFactor) =
    velocity(sensor(w), E, fc)

(w::Probe1d)(E::Real; kw...) = velocity(E; kw...) =
    velocity(sensor(w), E; kw...)

(w::Probe1d)(E::Real, fc::CorrFactor) =
    velocity(sensor(w), E, fc)

