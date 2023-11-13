-# Temperature/Pressure/etc correction

abstract type AbstractAnemCorrect end


resistance(c::AbstractAnemCorrect) = c.Rw
temperature(c::AbstractAnemCorrect) = c.Tw
reftemp(c::AbstractAnemCorrect) = c.T
flowfluid(c::AbstractAnemCorrect) = c.fluid


anemcorrectfactor(cal::AbstractAnemCorrect, Tw, Rw, T) =
    sqrt(resistance(cal)/Rw * (temperature(cal)-reftemp(cal)) / (Tw-T))


struct TempCorrect{U<:AbstractFloat,Fluid} <: AbstractAnemCorrect
    "Calibration operating resistance"
    Rw::U
    "Calibration operating resistance temperature"
    Tw::U
    "Calibration temperature"
    T::U
    "Calibration fluid"
    fluid::Fluid
    "Calibration kinematic viscosity"
    ν::U
end

anemcorrect(cal::TempCorrect, E, meas::TempCorrect) =
    E * anemcorrectfactor(cal, temperature(meas), resistance(meas), reftemp(meas))

anemcorrect(cal::TempCorrect, E, T, P, x, Rw, Tw) =
    E * anemcorrect(cal, Tw, Rw, T)

    
struct WireCorrect{U<:AbstractFloat,Fluid} <: AbstractAnemCorrect
    "Calibration operating resistance"
    Rw::U
    "Calibration operating resistance temperature"
    Tw::U
    "Calibration temperature"
    T::U
    "Calibration Fluid"
    fluid::Fluid
    "Calibration kinematic viscosity"
    ν::U
    "Calibration k⋅Prⁿ"
    ϕ::U
    "Prandtl number exponent"
    n::U
end

function anemcorrect(cal::WireCorrect, E, meas::WireCorrect)
    f = anemcorrectfactor(cal, temperature(meas), resistance(meas), reftemp(meas))
    return sqrt(cal.ϕ/meas.ϕ) * f
end


struct GlassbeadCorrect{U<:AbstractFloat,Fluid} <: AbstractAnemCorrect
    "Calibration operating resistance"
    Rw::U
    "Calibration operating resistance temperature"
    Tw::U
    "Calibration temperature"
    T::U
    "Calibration Fluid"
    fluid::Fluid
    "Calibration kinematic viscosity"
    ν::U
    "Calibration k⋅Prⁿ"
    ϕ::U
    "Prandtl number exponent"
    n::U
    "Heat conductivity coefficient of the shell"
    β::U
    "Glass bead A/L"
    c1::U
    "Leads N⋅√(γ/D⋅kAPf)"
    c2::U
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

function anemcorrect(corr::GlassbeadCorrect, E, Rw, Tw, ρ, μ, k, Pr, args...)
    
    Rwc = corr.Rw
    Twc = corr.Tw
    Tac = corr.Ta

    β = corr.β
    ϕc = corr.ϕ
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

