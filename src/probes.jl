
abstract type AbstractProbe end

abstract type AbstractProbe1d <: AbstractProbe end
abstract type AbstractProbe2d <: AbstractProbe end

abstract type AbstractProbe3d <: AbstractProbe end

struct HWSupport
    "Support model"
    model::String
    "Support tag (storage and control)"
    tag::String
    "Nominal resistance of the support"
    R::Float64
    "Diameter of the support"
    D::Float64
    "Length of the support"
    L::Float64
    "Length of support cables"
    Lc::Float64
end

"""
    `HWSupport(;model="55H21", tag="", R=0.44, D=4.0, L=235, Lc=765)`

Stores information about Probe support. Noty strictly necessary but 
could be useful for future reference

"""
HWSupport(;model="", tag="", R=0.0, D=0.0, L=0.0, Lc=0.0) = HWSupport(model, tag, R, D, L, Lc)

struct HWCable
    "Cable model/type"
    model::String
    "Cable tag (storage and control)"
    tag::String
    "Nominal resistance of the cable"
    R::Float64
    "Length of the cable"
    L::Float64
    "Connector types (left and right)"
    conn::{String,String}
end

"""
    `HWCable(model="", tag="", R=0.2, L=4.0, ("BNC","BNC"))`

Stores information about the cable used. Not strictly necessary but
could be useful.

"""
HWCable(;model="A1863 - 50Ω impedance", tag="", R=0.2, L=4.0, conn=("BNC", "BNC")) = HWCable(model, tag, R, L, conn)


struct HWBridge
    "Anemometer system model"
    model::String
    "Anemometer tag"
    tag::String
    "Specific bridge id"
    id::Int
    "Output offset"
    offset::Float64
    "Output gain"
    gain::Float64
    "Output low pass filter"
    lpfilt::Float64
    "Bridge ration"
    bridgeratio::String
    "Bridge amplifier gain"
    ampgain::Float64
    "Bridge filter"
    filter::Float64
    "Cable compensation"
    cablecomp::Int
    "Cooling interval"
    coolinterv=Float64
end

"""
    `HWBridge(model="", tag="", id=1, offset=0, gain=1, lpfilt=1e3,
              bridgeration="20:1", ampgain=8, filter=3, 
              cablecomp=0, coolinterv=2)`

Information about hotwire hardware system.

"""
HWBridge(;model="Streamline", tag="", id=1, offset=0.0, gain=1.0,
         lpfilt=1e3, bridgeratio="20:1", ampgain=8, filter=3,
         cablecomp=0, coolinterv=2.0) =
             HWBridge(model, tag, id, offset, gain, lpfilt,
                      bridgeratio, ampgain, filter, cablecomp,
                      coolinterv)

mutable struct Probe1d{Anem<:AbstractThermalAnemometer,Calibr} <: AbstractProbe1d
    "Probe model"
    model::String
    "Probe tag"
    tag::String
    "Probe leads resistance"
    RL::Float64
    "Anemometer sensor information"
    sensor::Anem
    "Calibration data"
    cal::Calibr
    "Hardware information"
    bridge::HWBridge
    "Probe support"
    support::HWSupport
    "Connecting cable"
    cable::HWCable
end

(anem::Probe1d)(E, temp) = anem.cal(anem.sensor, E, temp)
(anem::Probe1d)(E) = anem.cal(anem.sensor, E, anem.cal.T0)

mutable struct Probe2d{Anem<:AbstractThermalAnemometer,Calibr} <: AbstractProbe2d
    "Probe model"
    model::String
    "Probe tag"
    tag::String
    "Probe leads resistance"
    RL::NTuple{2,Float64}
    "Sensor information"
    sensor::NTuple{2,Anem}
    "Sensor calibration"
    cal::NTuple{2,Calibr}
    "k² directional calibration for each sensor"
    k²::NTuple{2,Float64}
    "h² directional calibration for each sensor"
    h²::NTuple{2,FLoat64}
    "Matrix with the cosine of each wire with respect to x, y and z axes"
    cosϕ::Matrix{Float64}
    "Hardware information"
    bridge::NTuple{2,HWBridge}
    "Probe support"
    support::HWSupport
    "Connecting cable"
    cable::NTuple{2,HWCable}
    
end


mutable struct Probe3d{Anem<:AbstractThermalAnemometer,Calibr} <: AbstractProbe3d
    "Probe model"
    model::String
    "Probe tag"
    tag::String
    "Probe leads resistance"
    RL::NTuple{3,Float64}
    "Sensor information"
    sensor::NTuple{3,Anem}
    "Sensor calibration"
    cal::NTuple{3,Calibr}
    "k² directional calibration for each sensor"
    k²::NTuple{3,Float64}
    "h² directional calibration for each sensor"
    h²::NTuple{3,FLoat64}
    "Matrix with the cosine of each wire with respect to x, y and z axes"
    cosϕ::Matrix{Float64}
    "Hardware information"
    bridge::NTuple{3,HWBridge}
    "Probe support"
    support::HWSupport
    "Connecting cable"
    cable::NTuple{3,HWCable}
    
end



    

