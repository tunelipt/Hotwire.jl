mutable struct Probe3d{Anem<:AbstractThermalAnemometer,Calibr,Setup} <: AbstractProbe3d
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
    "Anemometer setup"
    setup::Setup
    "Factor to calculate effective velocity from calibration velocity"
    c2e::SVector{3,Float64}
    "Matrix to calculate velocity in wire coordinates from effective vel."
    A::SMatrix{3,3,Float64,9}
end

"""
    `Probe3d(sensor, cal, k², h², sensor="")`

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
 - `setup`: hardware/system  configuration
 - `ang`: Angle that each wire makes with each coordinate axes
"""
function Probe3d(sensor, cal, k², h², setup=""; ang=[125.264 45 114.094;
                                                   125.264 135 114.094;
                                                   125.264 90 35.264])

    # Calculate the cosine of the angles
    cosϕ = SMatrix{3,3,Float64}(cosd.(ang))
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
