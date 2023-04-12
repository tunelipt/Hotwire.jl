mutable struct Probe3d{Anem<:AbstractThermalAnemometer,Calibr,Setup} <: AbstractProbe3d
    "Sensor information"
    sensor::NTuple{3,Anem}
    "Sensor calibration"
    cal::NTuple{3,Calibr}
    "k² directional calibration coefficient"
    k²::Vector{Float64}
    "h² directional calibration coefficient"
    h²::Vector{Float64}
    "Index of h² for each wire"
    ih::Vector{Int}
    "Cosine of wire direction with respect to axes x, y and z"
    cosϕ::Matrix{Float64}
    "Factor to calculate effective velocity from calibration velocity"
    c2e::Vector{Float64}
    "Matrix to calculate velocity in wire coordinates from effective vel."
    A::Matrix{Float64}
    "Index of axis aloing which calibration is carried out"
    idircal::Int
    "Anemometer setup"
    setup::Setup
end

const cosϕ_dantec =  [-sqrt(1/3)  1/sqrt(2)   -sqrt(1-1/3-1/2);
                      -sqrt(1/3)  -1/sqrt(2)  -sqrt(1-1/3-1/2);
                      -sqrt(1/3)      0.0         sqrt(1-1/3) ]

"""
`Probe3d(sensor, cal, k², h², setup=""; ih=[3,1,2], idircal=1,
         cosang=cosϕ_dantec)`

Creates a 3d probe. This completely characterizes a 3d probe. It includes
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

The function assumes that the calibration was carried along axis `idircal`.

### Arguments:

 - `sensor`: a tuple of `AbstractThermalAnemometer` objects (3)
 - `cal`: Tuple of `CalibrCurve` objects, one calibration for each sensor
 - `k²`: directional calibration coefficient k² for each wire. See [`dircalibr`](@ref).
 - `h²`: directional calibration coefficient k² for each wire. See [`dircalibr`](@ref).
 - `setup`: hardware/system  configuration
 - `ih`: index of velocity component where h² is applied for each wire
 - `idircal`: Index of axis along which the calibration is carried out
 - `cosang`: Cossine of angle that each wire makes with each coordinate axes

The angle that each wire makes with the probe coordinates system is a matrix
where each *row* `i` is the angle that wire `i` makes with axes x, y and z.

For Dantex triaxial probes, wires 1 and make an angle of 45° with axis y, wire
3 is in plane xz and cos²θ = 1/3 where θ is the angle that the wires 1 - 3 make
with coordinate axis x.
"""
function Probe3d(sensor, cal, k², h², setup=""; ih=[3,1,2], idircal=1,
                 cosang=cosϕ_dantec)
    K = ones(3,3)
    for i in 1:3
        K[i,i] = k²[i]
        K[i,ih[i]] = h²[i]
    end
    
    cosϕ = cosang
    c2e =  [(1.0 + k²[1] + h²[1])*cosϕ[1,idircal]^2,
            (1.0 + k²[2] + h²[2])*cosϕ[2,idircal]^2,
            (1.0 + k²[3] + h²[3])*cosϕ[3,idircal]^2 ]
    
    return Probe3d((sensor[1], sensor[2], sensor[3]),
                   (cal[1], cal[2], cal[3]), k², h², ih, cosϕ, 
                   c2e, inv(K), idircal, setup)
end

reftemp(anem::Probe3d) = reftemp(anem.sensor[1])

function velocity(anem::Probe3d, E₁, E₂, E₃, T)

    # Calibration curve and temperature correction for
    # each sensor and calculation of effective velocity for each wire
    Uc1 = anem.cal[1](anem.sensor[1], E₁, T)
    Uc2 = anem.cal[2](anem.sensor[2], E₂, T)
    Uc3 = anem.cal[3](anem.sensor[3], E₃, T)
    Uᵉ = [Uc1^2 * anem.c2e[1], 
          Uc2^2 * anem.c2e[2], 
          Uc3^2 * anem.c2e[3]]
    
    # Calculate speed in wire coordinates
    Uʷ = anem.A * Uᵉ
    
    U1 = sqrt(abs(Uʷ[1]))
    U2 = sqrt(abs(Uʷ[2]))
    U3 = sqrt(abs(Uʷ[3]))
    
    return (-U1*anem.cosϕ[1,1] - U2*anem.cosϕ[2,1] - U3*anem.cosϕ[3,1],
            -U1*anem.cosϕ[1,2] - U2*anem.cosϕ[2,2] - U3*anem.cosϕ[3,2],
            -U1*anem.cosϕ[1,3] - U2*anem.cosϕ[2,3] - U3*anem.cosϕ[3,3])
end

(anem::Probe3d)(E₁, E₂, E₃, T) = velocity(anem, E₁, E₂, E₃, T)
(anem::Probe3d)(E₁, E₂, E₃) = velocity(anem, E₁, E₂, E₃, reftemp(anem.sensor[1]))


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
    idir = anem.idircal
    # Efective cooling velocity from probe calibration
    u1cos = Uc1 .^ 2 * cosw[1,idir]^2
    u2cos = Uc2 .^ 2 * cosw[2,idir]^2
    u3cos = Uc3 .^ 2 * cosw[3,idir]^2
    
    Ak1 = U₁² .- u1cos;  Ah1 = U₃² .- u1cos; y1 = u1cos .- U₂²
    Ak2 = U₂² .- u2cos;  Ah2 = U₁² .- u2cos; y2 = u2cos .- U₃²
    Ak3 = U₃² .- u3cos;  Ah3 = U₂² .- u3cos; y3 = u3cos .- U₁²

    
    fit1 = linearfit(y1, Ak1, Ah1) 
    fit2 = linearfit(y2, Ak2, Ah2) 
    fit3 = linearfit(y3, Ak3, Ah3) 

    k = [fit1[1], fit2[1], fit3[1]]
    h = [fit1[2], fit2[2], fit3[2]]
    
    return k,h
    

end

