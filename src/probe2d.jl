

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


function dircalibr(anem::Probe2d, ang, Uc, E1, E2, T)
    Ux = Uc * cosd(ang)
    Uy = Uc * sind(ang)

    
end



    
    
    
