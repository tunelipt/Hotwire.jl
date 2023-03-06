
using StaticArrays

abstract type AbstractProbe end

abstract type AbstractProbe1d <: AbstractProbe end
abstract type AbstractProbe2d <: AbstractProbe end

abstract type AbstractProbe3d <: AbstractProbe end


mutable struct Probe1d{Anem<:AbstractThermalAnemometer,Calibr,Setup} <: AbstractProbe1d
    "Anemometer sensor information"
    sensor::Anem
    "Calibration data"
    cal::Calibr
    "Anemometer setup"
    setup::Setup
end

"""
    `Probe1d(sensor, cal, setup)`

Creates a 1d probe. This completely characterizes a 1d probe. It includes
information on calibration, bridge configuration, output filter, etc.

# Arguments:

 - `sensor`: an `AbstractThermalAnemometer` object
 - `cal`: `CalibrCurve` object

"""
function Probe1d(sensor, cal, setup="")
    return Probe1d(sensor, cal, setup)
end

                                                        
(anem::Probe1d)(E, temp) = anem.cal(anem.sensor, E, temp)
(anem::Probe1d)(E) = anem.cal(anem.sensor, E, anem.cal.T0)


"""
    `velocity!(anem, E, T, 1)`

Calculate, in place, the velocity from an anemometer. In this function, 
the input is an abstract array where 1 column with the bridge output is given by argument `idx` represent the output of the anemometer. A unique temperature should be used.

"""
function velocity!(anem::Probe1d, E::AbstractMatrix, T, idx=1)

    npts = size(E,1)
    for i in 1:npts
        U = anem(E[i,idx], T)
        E[i,idx] = U
    end
end

    
"""
    `velocity(anem, E, T, 1)`

Calculate the velocity from an anemometer. In this function, 
the input is an abstract array where 1 column with the bridge output is given by argument `idx` represent the output of anemometer. A unique temperature should be used

"""
function velocity(anem::Probe1d, E::AbstractMatrix, T, idx=1)

    npts = size(E,1)
    U = zeros(npts)
    for i in 1:npts
        U[i] = anem(E[i,idx], T)
    end
    return U
end

    



