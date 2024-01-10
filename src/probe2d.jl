import Statistics: mean


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

function vel_aux(w::Probe2d, Uc1::Real, Uc2::Real)

    k₁², k₂² = w.k²

    a = 1 / (2*(1 - k₁²*k₂²))
    
    U₁ = sqrt( a * (Uc2^2 * (1 + k₂²) - Uc1^2 * k₂² * (1 + k₁²) ) )
    U₂ = sqrt( a * (Uc1^2 * (1 + k₁²) - Uc2^2 * k₁² * (1 + k₂²) ) )

    Ux = 1/sqrt(2) * (U₁ + U₂)
    Uy = 1/sqrt(2) * (U₁ - U₂)

    return Ux, Uy

end

function velf(w::Probe2d, E1::Real, E2::Real)
    Uc1 = velf(sensor(w,1),E1)
    Uc2 = velf(sensor(w,2),E2)

    return vel_aux(w, Uc1, Uc2)
end

function velocity(w::Probe2d, E1::Real, E2::Real, fc::CorrFactor)
    w1 = sensor(w,1); w2 = sensor(w,2)
    E1c = outvolt(w1, sensorvolt(w1, E1)*fc.f[1])
    E2c = outvolt(w2, sensorvolt(w2, E2)*fc.f[2])
    
    rnu = fc.nu[1]/kinvisc(w1) # We will assume the others are the same...
    Ux,Uy = velf(w, E1c, E2c)
    return Ux * rnu, Uy*rnu
end

velocity(w::Probe2d, E::NTuple, fc::CorrFactor) =
    velocity(w, E[1], E[2], fc)

velocity(w::Probe2d, E1::Real, E2::Real; kw...) = 
    velocity(w, E1, E2, correct(w, E1, E2; kw...))

velocity(w::Probe2d, E::NTuple; kw...) = 
    velocity(w, E[1], E[2]; kw...)

function velocity!(U::AbstractMatrix{T}, w::Probe2d,
                   E::AbstractMatrix{T}, fc::CorrFactor) where {T}
    @assert size(U,1) == size(E,1)
    @assert size(U,2) == size(E,2) == 2

    i = firstindex(U)
    for ee in eachrow(E)
        ux,uy = velocity(w, ee[1], ee[2], fc)
        U[i,1] = ux
        U[i,2] = uy
        i = i + 1
    end
    return U
end
function velocity!(U::AbstractMatrix{T}, w::Probe2d,
                   E::AbstractMatrix{T}; kw...) where {T}
    Em = mean(E, dims=1)
    return velocity!(U, w, E, correct(w, Em[1,1], Em[1,2], Em[1,3]; kw...))
end

velocity(w::Probe2d, E::AbstractMatrix, fc::CorrFactor) =
    velocity!(zeros(size(x)...), w, E, fc)

velocity(w::Probe2d, E::AbstractMatrix; kw...) =
    velocity!(zeros(size(x)...), w, E; kw...)

(w::Probe2d)(E1, E2, fc::CorrFactor) = velocity(w, E1, E2, fc)
(w::Probe2d)(E1, E2; kw...) = velocity(w, E1, E2; kw...)
(w::Probe2d)(E::NTuple, fc::CorrFactor) = velocity(w, E[1], E[2], fc)
(w::Probe2d)(E::NTuple; kw...) = velocity(w, E[1], E[2]; kw...)
(w::Probe2d)(E::AbstractMatrix, fc::CorrFactor) = velocity(w, E, fc)
(w::Probe2d)(E::AbstractMatrix; kw...) = velocity(w, E; kw...)



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



    
    
    
