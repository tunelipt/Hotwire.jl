# Temperature/Pressure/etc correction

abstract type AbstractAnemometerCorrection end

abstract type AbstractCTACorrection <: AbstractAnemometerCorrection
abstract type AbstractCCACorrection <: AbstractAnemometerCorrection


struct TempCorrection{T,U,V} <: AbstractCTACorrection
    Rw::T
    Tw::U
    Ta::V
end


function correct(cta::CTASensor, E; T=C.Ta, R=C.Rw, ...) where
    if R != C.Rw
        Tw1 = temperature(Rsens, R)
    else
        Tw1 = C.Tw
    end
    
    return sqrt(C.Rw/R * (R.tw - R.Ta) / (Tw1 - T))
end

struct WireCorrection{T,U,V,W,X,Y,Z} <: AbstractCTACorrection
    Rw::T
    Tw::U
    Ta::V
    ρ::W
    μ::X
    k::Y
    Pr::Z
    n
end

function correct(Rsens::AbstractResistor, C::WireCorrection;
                 R=C.Rw, T=C.Ta, rho=C.ρ, mu=C.μ, k=C.k, Pr=C.Pr)
                 

