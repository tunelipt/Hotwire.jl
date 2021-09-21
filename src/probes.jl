
using StaticArrays

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
    conn::Tuple{String,String}
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
    coolinterv::Float64
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
    "Anemometer sensor information"
    sensor::Anem
    "Calibration data"
    cal::Calibr
    "Probe model"
    model::String
    "Probe tag"
    tag::String
    "Probe leads resistance"
    RL::Float64
    "Hardware information"
    bridge::HWBridge
    "Probe support"
    support::HWSupport
    "Connecting cable"
    cable::HWCable
end

"""
    `Probe1d(sensor, cal; model="model", tag="HW 1d", RL=0.5, 
             bridge=HWBridge(), support=HWSupport(), cable=HWCable())`

Creates a 1d probe. This completely characterizes a 1d probe. It includes
information on calibration, bridge configuration, output filter, etc.

# Arguments:

 - `sensor`: an `AbstractThermalAnemometer` object
 - `cal`: `CalibrCurve` object
 - `model`: Probe model
 - `tag`: Probe id tag
 - `RL`: Lead resistance
 - `bridge`: Bridge configuration
 - `support`: `HWSupport` object specifying the probe support
 - `cable`: `HWCable` object with information on cable

"""
function Probe1d(sensor, cal; model="", tag="", RL=0.5,
                 bridge=HWBridge(id=1),
                 support=HWSupport(), cable=HWCable())
    return Probe1d(sensor, cal, model, tag, RL, bridge, support, cable)
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
        U = anem(E[i,idx], T])
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
        U[i] = anem(E[i,idx], T])
    end
    return U
end

    

mutable struct Probe2d{Anem<:AbstractThermalAnemometer,Calibr} <: AbstractProbe2d
    "Sensor information"
    sensor::NTuple{2,Anem}
    "Sensor calibration"
    cal::NTuple{2,Calibr}
    "k² directional calibration for each sensor"
    k²::NTuple{2,Float64}
    "h² directional calibration for each sensor"
    h²::NTuple{2,Float64}
    "Matrix with the cosine of each wire with respect to x, y and z axes"
    cosϕ::Matrix{Float64}
    "Probe model"
    model::String
    "Probe tag"
    tag::String
    "Probe leads resistance"
    RL::NTuple{2,Float64}
    "Hardware information"
    bridge::NTuple{2,HWBridge}
    "Probe support"
    support::HWSupport
    "Connecting cable"
    cable::NTuple{2,HWCable}
    
end

mutable struct Probe3d{Anem<:AbstractThermalAnemometer,Calibr} <: AbstractProbe3d
    "Sensor information"
    sensor::NTuple{3,Anem}
    "Sensor calibration"
    cal::NTuple{3,Calibr}
    "k² directional calibration for each sensor"
    k²::NTuple{3,Float64}
    "h² directional calibration for each sensor"
    h²::NTuple{3,Float64}
    "Cosine of wire direction with respect to axes x, y and z"
    cosϕ::SMatrix{3,3,Float64,9}
    "Probe model"
    model::String
    "Probe tag"
    tag::String
    "Probe leads resistance"
    RL::NTuple{3,Float64}
    "Factor to calculate effective velocity from calibration velocity"
    c2e::SVector{3,Float64}
    "Matrix to calculate velocity in wire coordinates from effective vel."
    A::SMatrix{3,3,Float64,9}
    "Hardware information"
    bridge::NTuple{3,HWBridge}
    "Probe support"
    support::HWSupport
    "Connecting cable"
    cable::NTuple{3,HWCable}
end

"""
    `Probe3d(sensor, cal, k², h²; ϕ=ϕ, model="model", tag="HW 1d", RL=0.5, 
             bridge=br, support=su, cable=cc)`

Creates a 2d probe. This completely characterizes a 1d probe. It includes
information on calibration, directional calibration, bridge configuration,
 output filter, etc.
The effective velocity acting on a wire is given by 
``
(Uᵉ)² = U₁² + h²⋅U₂² + k²⋅U₃²
``
In this equation, Uᵉ is the effective velocity, U₁ is the velocity normal 
to the wire and in the direction of the probe body, U₂ is the velocity 
normal to the flow but also normal to the body and U₃ is the velocity 
along the wire.

The function assumes that the calibration was carried out with flow along
the x axis. 

# Arguments:

 - `sensor`: a tuple of `AbstractThermalAnemometer` objects (3)
 - `cal`: `CalibrCurve` object
 - `model`: Probe model
 - `tag`: Probe id tag
 - `RL`: Lead resistance
 - `bridge`: Bridge configuration
 - `support`: `HWSupport` object specifying the probe support
 - `cable`: `HWCable` object with information on cable

"""
function Probe3d(sensor, cal, k², h²; ϕ=[125.264 45 114.094;
                                         125.264 135 114.094;
                                         125.264 90 35.264],
                 model="55R91", tag="", RL=0.5,
                 bridge=(HWBridge(), HWBridge(), HWBridge()),
                 support=HWSupport(),
                 cable=(HWCable(), HWCable(), HWCable()))

    # Calculate the cosine of the angles
    cosϕ = SMatrix{3,3,Float64}(cosd.(ϕ))
    c2e = SVector{3,Float64}( (1.0 + k²[1] + h²[1])*cosϕ[1,1]^2,
                              (1.0 + k²[2] + h²[2])*cosϕ[2,1]^2,
                              (1.0 + k²[3] + h²[3])*cosϕ[3,1]^2 )
    A = SMatrix{3,3,Float64}([k²[1] 1.0      h²[1];
                              h²[2] k²[2]    1.0;
                              1.0   h²[3]  k²[3]])
    return Probe3d((sensor[1], sensor[2], sensor[3]),
                   (cal[1], cal[2], cal[3]), (k²[1], k²[2], k²[3]),
                   (h²[1], h²[2], h²[3]), cosϕ, model, tag, (RL,RL,RL),
                   c2e, inv(A), (bridge[1],bridge[2],bridge[3]),
                   support, (cable[1],cable[2],cable[3]))
end

function (anem::Probe3d)(E₁, E₂, E₃, T)

    # Calibration curve and temperature correction for
    # each sensor and calculation of effective velocity for each wire
    Uᵉ = SVector(anem.cal[1](anem.sensor[1], E₁, T)^2 * anem.c2e[1], 
                 anem.cal[2](anem.sensor[2], E₂, T)^2 * anem.c2e[2], 
                 anem.cal[3](anem.sensor[3], E₃, T)^2 * anem.c2e[3])
    
    # Calculate speed in wire coordinates
    Uʷ = anem.A * Uᵉ
    
    U1 = sqrt(max(Uʷ[1], 0.0))
    U2 = sqrt(max(Uʷ[2], 0.0))
    U3 = sqrt(max(Uʷ[3], 0.0))
    
    return (-U1*anem.cosϕ[1,1] - U2*anem.cosϕ[2,1] - U3*anem.cosϕ[3,1],
            -U1*anem.cosϕ[1,2] - U2*anem.cosϕ[2,2] - U3*anem.cosϕ[3,2],
            -U1*anem.cosϕ[1,3] - U2*anem.cosϕ[2,3] - U3*anem.cosϕ[3,3])
end


"""
    `velocity!(anem, E, T, 1:3)`

Calculate, in place, the velocity from an anemometer. In this function, 
the input is an abstract array where 3 columns given by argument `idx` 
represent the output of each anemometer. A unique temperature should be 
used

"""
function velocity!(anem::Probe3d, E::AbstractMatrix, T, idx=1:3)

    npts = size(E,1)
    for i in 1:npts
        U, V, W = anem(E[i,idx[1]], E[i,idx[2]], E[i,idx[3]], T)
        E[i,idx[1]], E[i,idx[2]], E[i,idx[3]] = U, V, W
    end

end

"""
    `velocity(anem, E, T, 1:3)`

Velocity calculation for arrays. A new array where each column represents
a component of the velocity is created

"""
function velocity(anem::Probe3d, E::AbstractMatrix, T, idx=1:3)
    npts = size(E,1)
    
    Uout = zeros(npts, 3)
    
    for i in 1:npts
        U, V, W = anem(E[i,idx[1]], E[i,idx[2]], E[i,idx[3]])
        Uout[i,idx[1]], Uout[i,idx[2]], Uout[i,idx[3]] = U, V, W
    end

    return Uout
end

