

mutable struct Probe2d{Anems,Setup} <: AbstractProbe2d
    "Sensor information"
    sensor::Anems
    "k² directional calibration for each sensor"
    k²::NTuple{2,Float64}
    "Anemometer setup"
    setup::Setup
end

sensor(w::Probe2d) = w.sensor
sensor(w::Probe2d, i::Integer) = w.sensor[i]

reftemp(anem::Probe2d) = reftemp(anem.sensor[1])

function correct(w::Probe2d, E1::Real, E2::Real; kw...)
    fc1 = correct(sensor(w,1), E1, kw...)
    fc2 = correct(sensor(w,2), E2, kw...)
    return CorrFactor((fc1.f, fc2.f), (fc1.nu, fc2.nu))
end

function velocity(w::Probe2d, E1::Real, E2::Real, fc::CorrFactor)
    
    
    Uc1 = sensor(w,1).fit(E1*fc.f[1])
    Uc2 = sensor(w,2).fit(E2*fc.f[2])
    

    k₁², k₂² = w.k²

    a = 1 / (2*(1 - k₁²*k₂²))
    
    U₁ = sqrt( a * (Uc2^2 * (1 + k₂²) - Uc1^2 * k₂² * (1 + k₁²) ) )
    U₂ = sqrt( a * (Uc1^2 * (1 + k₁²) - Uc2^2 * k₁² * (1 + k₂²) ) )

    Ux = 1/sqrt(2) * (U₁ + U₂)
    Uy = 1/sqrt(2) * (U₁ - U₂)

    ν_cal = kinvisc(sensor(w,1))
    fre = fc.nu[1] / ν_cal
    return Ux*fre, Uy*fre
end

function velocity(w::Probe2d, E1::Real, E2::Real; kw...)
    fc = correct(w, E1, E2; kw...)
    return velocity(w, E1, E2, fc)
end

(w::Probe2d)(E1::Real, E2::Real, fc::CorrFactor) = velocity(w, E1, E2, fc)
(w::Probe2d)(E1::Real, E2::Real; kw...) = velocity(w, E1, E2; kw...)



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



    
    
    
