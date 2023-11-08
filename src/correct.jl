# Temperature/Pressure/etc correction

abstract type AbstractAnemCorrect end


anemcorrectfactor(Tac, Twc, Rwc, Ta, Tw, Rw) = sqrt(Rwc/Rw * (Twc-Ta) / (Tw-Ta))


struct TempCorrect{T<:AbstractFloat} <: AbstractAnemCorrect
    "Calibration operating resistance"
    Rw::T
    "Calibration operating resistance temperature"
    Tw::T
    "Calibration temperature"
    Ta::T
    "Calibration kinematic viscosity"
    ν::T
end



function anemcorrect(corr::TempCorrect, Rsens::AbstractResistor, E,
                     T, P=0, R=corr.Rw, x=0)
    Rwc = corr.Rw
    if R != Rwc
        Rw = R
        Tw = temperature(Rsens, R)
    else
        Rw = Rwc
        Tw = corr.Tw
    end

    return E * anemcorrectfactor(corr.Ta, corr.Tw, Rwc, T, Tw, Rw)
end

function anemcorrect(corr::TempCorrect, Rsens::AbstractResistor, E;
                     T=corr.Ta, R=corr.Rw, kw...)
    return E * anemcorrectfactor(corr, Rsens, T, 0, R, 0)
end

    
struct WireCorrect{T<:AbstractFloat} <: AbstractAnemCorrect
    "Calibration operating resistance"
    Rw::T
    "Calibration operating resistance temperature"
    Tw::T
    "Calibration temperature"
    Ta::T
    "Calibration pressure"
    Pa::T
    "Calibration kinematic viscosity"
    ν::T
    "Calibration k⋅Prⁿ"
    kPrn::T
    "Prandtl number exponent"
    n::T
end

function anemcorrect(corr::WireCorrect, Rsens::AbstractResistor, E, 
                     T=corr.Ta, P=corr.Pa, R=corr.Rw, x=corr.fluid)
    Rwc = corr.Rw
    Twc = corr.Tw
    Tac = corr.Ta

    if R != Rwc
        Rw = R
        Tw = temperature(Rsens, R)
    else
        Rw = Rwc
        Tw = Twc
    end

    ρ, μ, k, Pr = fluidprops(fluid, T, P)

    kPrn = k * Pr^corr.n

    return sqrt(corr.kPrn/kPrn) * anemcorrectfactor(Tac, Twc, Rwc, T, Tw, Rw)
end


struct GlassbeadCorrect{T<:AbstractFloat} <: AbstractAnemCorrect
    "Calibration operating resistance"
    Rw::T
    "Calibration operating resistance temperature"
    Tw::T
    "Calibration temperature"
    Ta::T
    "Calibration pressure"
    Pa::T
    "Calibration kinematic viscosity"
    ν::T
    "Calibration k⋅Prⁿ"
    ϕ::T
    "Prandtl number exponent"
    n::T
    "Heat conductivity coefficient of the shell"
    β::T
    "Glass bead A/L"
    c1::T
    "Leads N⋅√(γ/D⋅kAPf)"
    c2::T
end

function mf58correct(Rsens ;R=1e2, T=298.15, P=93e3, n=1/3, N=2, beta=0.2, fluid=Air)
    Tw = temperature(Rsens, R)
    Rw = R

    ρ, μ, ka, Pr = fluidprops(fluid, T, P)
    ν = μ/ρ
    
    ϕ = Pr^n * ka
    D = 2.0e-3
    L = 4.0e-3
    A = π*D*L
    c1 = A / D
    d = 0.5e-3
    P = π*d
    Af = π*d^2/4
    kf = 30.0 #W/mK
    γ = sqrt(d/D)
    β = 150.0 # W/K
    c2 = N * sqrt(γ*kf*Af*P/D)
    return GlassbeadCorrect(Rw, Tw, T, P, ν, ϕ, n, β, c1, c2)
end   
    
function anemcorrect(corr::GlassbeadCorrect, Rsens::AbstractResistor, E, 
                     T, P=corr.Pa, R=corr.Rw, x=Air)
    
    Rwc = corr.Rw
    Twc = corr.Tw
    Tac = corr.Ta
    if R != Rwc
        Rw = R
        Tw = temperature(Rsens, R)
    else
        Rw = Rwc
        Tw = Twc
    end
    β = corr.β
    ϕc = corr.ϕ
    ρ, μ, k, Pr = fluidprops(x, T, P)
    ϕ = k * Pr^corr.n
    
    Y = E*E / (Rw * (Tw - T))
    X = Y / (1-β*Y)

    c1 = corr.c1
    c2 = corr.c2

    u = (-c2 + sqrt(c2*c2 + 4*c1*X)) / 2c1

    fc = u*u / ϕ
    
    Xc = c1*fc*ϕc + c2*sqrt(fc*ϕc)
    Yc = Xc / (1 + β*Xc)
    
    return sqrt(Yc * Rwc * (Twc - Tac))
    
end

anemcorrect(corr::GlassbeadCorrect, Rsens::AbstractResistor, E;
            T=corr.Ta, P=corr.Pa, R=corr.Rw, x=Air) = anemcorrect(corr, Rsens, E,
                                                                  T, P, R, x)
