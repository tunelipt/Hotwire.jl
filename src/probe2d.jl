import Statistics: mean


mutable struct Probe2d{Anems,Setup,Device} <: AbstractProbe2d
    "Sensor information"
    sensor::Anems
    "k² directional calibration for each sensor"
    k²::NTuple{2,Float64}
    "Anemometer setup"
    setup::Setup
    dev::Device
end

sensor(w::Probe2d) = w.sensor
sensor(w::Probe2d, i::Integer) = w.sensor[i]




function vel_aux(w::Probe2d, Uc1::Real, Uc2::Real)

    k₁², k₂² = w.k²

    a = 1 / (2*(1 - k₁²*k₂²))
    
    U₁ = sqrt( a * (Uc2^2 * (1 + k₂²) - Uc1^2 * k₂² * (1 + k₁²) ) )
    U₂ = sqrt( a * (Uc1^2 * (1 + k₁²) - Uc2^2 * k₁² * (1 + k₂²) ) )

    Ux = 1/sqrt(2) * (U₁ + U₂)
    Uy = 1/sqrt(2) * (U₁ - U₂)

    return Ux, Uy

end

function vel_aux!(w::Probe2d,
                  Ux::AbstractArray, Uy::AbstractArray)

    k₁², k₂² = w.k²
    a = 1 / (2*(1 - k₁²*k₂²))
    for (i,Uc1,Uc2) in zip(eachindex(Uc1), Ux, Uy)
        U₁ = sqrt( a * (Uc2^2 * (1 + k₂²) - Uc1^2 * k₂² * (1 + k₁²) ) )
        U₂ = sqrt( a * (Uc1^2 * (1 + k₁²) - Uc2^2 * k₁² * (1 + k₂²) ) )
        
        Ux[i] = 1/sqrt(2) * (U₁ + U₂)
        Uy[i] = 1/sqrt(2) * (U₁ - U₂)
        
    end
    return nothing
    
end


function velocity(w::Probe2d, E1::Real, E2::Real; kw...)
    w1 = sensor(w,1); w2 = sensor(w,2)
    Uc1 = velocity(w1, E1; kw...)
    Uc2 = velocity(w2, E2; kw...)

    return vel_aux(w, Uc1, Uc2)
end

function velocity!(Ux::AbstractArray, Uy::AbstractArray,
                   w::Probe2d, E1::AbstractArray, E2::AbstractArray;
                   kw...)
    @assert size(E1) == size(E2) == size(Ux) == size(Uy)
    w1 = sensor(w,1); w2 = sensor(w,2)
    velocity!(Ux, w1, E1; kw...)
    velocity!(Uy, w2, E2; kw...)
    return Ux, Uy 
end

    

(probe::Probe2d)(E1, E2; kw...) = velocity(probe, E1, E2; kw...)




import LinearAlgebra: dot

"""
`dircalibr(ang, Uc, Uc1, Uc2)`

Perform directional calibration of a 2d probe (X wire).

### Arguments

 * `anem`: a 2d probe object
 * `ang`: Directional calibration angle in degrees
 * `Uc`: Directional calibration *constant* velocity
 * `Uc1`: probe velocity calibration output for wire 1
 * `Uc2`: probe velocity calibration output for wire 2

### Returns

Coefficients `k₁²` and `k₂²` as a tuple.
"""
function dircalibr(ang, Uc, Uc1, Uc2)
    
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



    
    
    
