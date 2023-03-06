

mutable struct Probe2d{Anem<:AbstractThermalAnemometer,Calibr,Setup} <: AbstractProbe2d
    "Sensor information"
    sensor::NTuple{2,Anem}
    "Sensor calibration"
    cal::NTuple{2,Calibr}
    "k² directional calibration for each sensor"
    k²::NTuple{2,Float64}
    "Anemometer setup"
    setup::Setup
end

function velocity(anem::Probe2d, E1, E2, T)

    Uc1 = anem.cal[1](E1, T)  # Calibration velocity of wire 1
    Uc2 = anem.cal[2](E2, T)  # Calibration velocity of wire 1

    k₁², k₂² = anem.k²

    a = 1 / (2*(1 - k₁²*k₂²))
    
    U₁ = sqrt( a * (Uc2^2 * (1 + k₂²) - Uc1^2 * k₂² * (1 + k₁²) ) )
    U₂ = sqrt( a * (Uc1^2 * (1 + k₁²) - Uc2^2 * k₁² * (1 + k₂²) ) )

    Ux = 1/sqrt(2) * (U₁ + U₂)
    Uy = 1/sqrt(2) * (U₁ - U₂)

    return Ux, Uy
end


import LinearAlgebra: dot

function dircalibr(anem::Probe2d, ang, Uc, Uc1, Uc2)
    
    A = @. (cosd(ang) + sind(ang))^2
    B = @. (cosd(ang) - sind(ang))^2

    X2A = @. (Uc^2 * A - Uc2^2)
    X2B = @. (Uc^2 * B - Uc2^2) #-

    X1A = @. (Uc^2 * A - Uc1^2)
    X1B = @. (Uc^2 * B - Uc1^2)

    # Fit the curves
    # k1*X1A = -X1B
    # k2*X2B = -X2A

    k1 = -dot(X1B, X1A) / dot(X1A, X1A)
    k2 = -dot(X2B, X2A) / dot(X2B, X2B)

    return k1, k2
end



    
    
    
