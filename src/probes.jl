
abstract type AbstractProbe end

abstract type AbstractProbe1d <: AbstractProbe end
abstract type AbstractProbe2d <: AbstractProbe end

abstract type AbstractProbe3d <: AbstractProbe end


struct Probe1d{Anem<:AbstractThermalAnemometer,Setup} <: AbstractProbe1d
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

resistor(w::Probe1d) = resistor(w.sensor)
sensor(w::Probe1d) = w.sensor

gain(w::Probe1d) = gain(w.sensor.signal)
offset(w::Probe1d) = offset(w.sensor.signal)
sensorvolt(w::Probe1d, E) = sensorvolt(w.sensor, E)
outsignal(w::Probe1d, E) = outsignal(w.sensor, E)

reftemp(w::Probe1d) = reftemp(w.sensor)
caltemp(w::Probe1d) = caltemp(w.sensor)
resistance(w::Probe1d) = resistance(w.sensor)
pressure(w::Probe1d) = pressure(w.sensor)
fluid(w::Probe1d) = fluid(w.sensor)

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

correct(w::Probe1d, E::Real; kw...) = correct(sensor(w), E; kw...)

velf(w::Probe1d, E) = velf(sensor(w),E)


velocity(w::Probe1d, E; kw...) = velocity(sensor(w), E; kw...)
velocity(w::Probe1d, E, fc::CorrFactor) =
    velocity(sensor(w), E, fc)

velocity!(U::AbstractVector, w::Probe1d, E::AbstractVector; kw...) =
    velocity!(U, sensor(w), E; kw...)
velocity(w::Probe1d, E::AbstractVector; kw...) = velocity(sensor(w), E; kw...)


(w::Probe1d)(E; kw...) = velocity(sensor(w), E; kw...)
(w::Probe1d)(E, fc::CorrFactor) =
    velocity(sensor(w), E, fc)


