mutable struct Probe3d{Anem<:AbstractThermalAnemometer,Calibr,Setup} <: AbstractProbe3d
    "Sensor information"
    sensor::NTuple{3,Anem}
    "Sensor calibration"
    cal::NTuple{3,Calibr}
    "Matrix with h² and k² for each wire"
    kh::SMatrix{3,3,Float64,9}
    "Cosine of wire direction with respect to axes x, y and z"
    cosϕ::SMatrix{3,3,Float64,9}
    "Factor to calculate effective velocity from calibration velocity"
    c2e::SVector{3,Float64}
    "Matrix to calculate velocity in wire coordinates from effective vel."
    A::SMatrix{3,3,Float64,9}
    "Anemometer setup"
    setup::Setup
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

The angle that each wire makes with the probe coordinates system is a matrix
where each *row* `i` is the angle that wire `i` makes with axes x, y and z.

For Dantex triaxial probes, wires 1 and make an angle of 45° with axis y, wire
3 is in plane xz and cos²θ = 1/3 where θ is the angle that the wires 1 - 3 make
with coordinate axis x.
"""
function Probe3d(sensor, cal, k², h², setup=""; ang=[125.2643897 45 114.0953355;
                                                     125.2643897 135 114.0953355;
                                                     125.2643897 90 35.2643897])
    # ang: angle 
    #   ang=[125.264 45 114.094;  125.264 135 114.094;   125.264 90 35.264])
    # DANTEC probes
    kh = @SMatrix [k²[1]   1.0   h²[1]
                   h²[2]  k²[2]    1.0
                    1.0   h²[3]  k²[3]]
    
    # Calculate the cosine of the angles
    cosϕ = SMatrix{3,3,Float64}(cosd.(ang))
    c2e = @SVector  [(1.0 + k²[1] + h²[1])*cosϕ[1,1]^2,
                     (1.0 + k²[2] + h²[2])*cosϕ[2,1]^2,
                     (1.0 + k²[3] + h²[3])*cosϕ[3,1]^2 ]
    
    return Probe3d((sensor[1], sensor[2], sensor[3]),
                   (cal[1], cal[2], cal[3]), kh, cosϕ, 
                   c2e, inv(kh), setup)
end

reftemp(anem::Probe3d) = reftemp(anem.sensor[1])

function velocity(anem::Probe3d, E₁, E₂, E₃, T)

    # Calibration curve and temperature correction for
    # each sensor and calculation of effective velocity for each wire
    Uᵉ = @SVector [anem.cal[1](anem.sensor[1], E₁, T)^2 * anem.c2e[1], 
                   anem.cal[2](anem.sensor[2], E₂, T)^2 * anem.c2e[2], 
                   anem.cal[3](anem.sensor[3], E₃, T)^2 * anem.c2e[3]]
    
    # Calculate speed in wire coordinates
    Uʷ = anem.A * Uᵉ
    
    U1 = sqrt(max(Uʷ[1], 0.0))
    U2 = sqrt(max(Uʷ[2], 0.0))
    U3 = sqrt(max(Uʷ[3], 0.0))
    
    return (-U1*anem.cosϕ[1,1] - U2*anem.cosϕ[2,1] - U3*anem.cosϕ[3,1],
            -U1*anem.cosϕ[1,2] - U2*anem.cosϕ[2,2] - U3*anem.cosϕ[3,2],
            -U1*anem.cosϕ[1,3] - U2*anem.cosϕ[2,3] - U3*anem.cosϕ[3,3])
end

(anem::Probe3d)(E₁, E₂, E₃, T) = velocity(anem, E₁, E₂, E₃, T)
(anem::Probe3d)(E₁, E₂, E₃) = velocity(anem, E₁, E₂, E₃, reftemp(anem.sensor[1]))


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


"""
`dircalibr(anem, Uc, Uc1, Uc2, Uc3, ang, angtilt)`

Performs an angular calibration of a triaxial probe.

### Arguments

* `anem`: `Probe3d` objetc
* `Uc`: Angular calibration velocity
* `Uc1`: calibration velocity for wire 1
* `Uc2`: calibration velocity for wire 1
* `Uc3`: calibration velocity for wire 1
* `ϕ`: roll angles used for directional calibration in degrees
* `θ`: Inclination angle used for directional calibration in degrees

The calibration velocity for each wire is the velocity obtained from
the calibration curve already corrected for temperature changes.
"""
function dircalibr(anem::Probe3d, Uc, Uc1, Uc2, Uc3, ϕ, θ=30.0)

    Ux = Uc .* cosd.(θ) .+ 0.0 .* ϕ
    uperp = Uc .* sind.(θ)
    Uz = uperp .* sind.(ϕ)
    Uy = - uperp .* cosd.(ϕ)

    # Cossine of angles of wire with probe coordinate system
    cosw = anem.cosϕ
    
    # Project this in wire coordinates
    U₁² = (Ux * cosw[1,1] + Uy * cosw[1,2] + Uz * cosw[1,3]).^2
    U₂² = (Ux * cosw[2,1] + Uy * cosw[2,2] + Uz * cosw[2,3]).^2
    U₃² = (Ux * cosw[3,1] + Uy * cosw[3,2] + Uz * cosw[3,3]).^2

    # Efective cooling velocity from probe calibration
    u1cos = Uc1 .^ 2 * cosw[1,1]^2
    u2cos = Uc2 .^ 2 * cosw[2,1]^2
    u3cos = Uc3 .^ 2 * cosw[3,1]^2
    
    Ak1 = U₁² .- u1cos;  Ah1 = U₃² .- u1cos; y1 = u1cos .- U₂²
    Ak2 = U₂² .- u2cos;  Ah2 = U₁² .- u2cos; y2 = u2cos .- U₃²
    Ak3 = U₃² .- u3cos;  Ah3 = U₂² .- u3cos; y3 = u3cos .- U₁²

    fit1 = CurveFit.fit_linear_model([Ak1 Ah1], y1)
    fit2 = CurveFit.fit_linear_model([Ak2 Ah2], y2)
    fit3 = CurveFit.fit_linear_model([Ak3 Ah3], y3)

    k = [fit1[1], fit2[1], fit3[1]]
    h = [fit1[2], fit2[2], fit3[2]]
    
    return k,h
    

end

