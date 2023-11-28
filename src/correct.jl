# Temperature/Pressure/etc correction

export AbstractAnemCorrect, TempCorrect, WireCorrect, GlassbeadCorrect
export mf58correct
    




struct TempCorrect{U<:Real}
    "Calibration operating resistance"
    Rw::U
    "Calibration operating resistance temperature"
    Tw::U
    "Calibration temperature"
    T::U
end

resistance(c::TempCorrect) = c.Rw
temperature(c::TempCorrect) = c.Tw
reftemp(c::TempCorrect) = c.T

function TempCorrect(c::TempCorrect;
                     T=reftemp(c), Rw=resistance(c), Tw=temperature(c))
    TempCorrect(Rw, Tw, T)
end
anemcorrectfactor(cal::TemCorrect, op::TempCorrect) =
    sqrt(resistance(cal)/resistance(op) *
    (temperature(cal) - reftemp(cal)) / (temperature(op) - reftemp(op)))

anemcorrectfactor(cal::TempCorrect, Tw, Rw, T) =
    sqrt(resistance(cal)/Rw * (temperature(cal)-reftemp(cal)) / (Tw-T))


resistance(c::AbstractAnemCorrect) = c.Rw
temperature(c::AbstractAnemCorrect) = c.Tw
reftemp(c::AbstractAnemCorrect) = c.T
kinvisc(c::AbstractAnemCorrect) = c.ν

"""
    `AbstractAnemCorrect`

Abstract class that corrects the output of thermal anemometers due to changes
in fluid properties and electronic configurations

The output of a thermal anemometer is directly related to the heat transfer from a 
"""
abstract type AbstractAnemCorrect end


anemcorrect(E, cal::TempCorrect, meas::TempCorrect) =
    E*anemcorrectfactor(cal, temperature(meas), resistance(meas), reftemp(meas))

anemcorrect(E, cal::TempCorrect, T, P, x, Rw, Tw) = E*anemcorrect(cal, Tw, Rw, T)

    
struct WireCorrect{U<:Real,Fluid} <: AbstractAnemCorrect
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

function WireCorrect(c::WireCorrect;
                     T=reftemp(c), Rw=resistance(c), Tw=temperature(c),
                     P=101325.0, x=fluid(c))
    ρ = density(x, T, P)
    μ = viscosity(x, T, P)
    cₚ = specheat(x, T, P)
    k = heatcond(x, T, P)

    Pr = cₚ * μ / k
    ν  = μ / ρ
    ϕ = k * Pr^c.n
    
    WireCorrect(Rw, Tw, T, x, ν, ϕ, n)
end

function anemcorrect(E,cal::WireCorrect, meas::WireCorrect)
    f = anemcorrectfactor(cal, temperature(meas), resistance(meas), reftemp(meas))
    return E*sqrt(cal.ϕ/meas.ϕ) * f
end

function anemcorrect(E,cal::WireCorrect, T, P, x, Rw, Tw)
    ρ = density(x, T, P)
    μ = viscosity(x, T, P)
    cₚ = specheat(x, T, P)
    k = heatcond(x, T, P)

    Pr = cₚ * μ / k
    ν  = μ / ρ
    return anemcorrect(E, cal, WireCorrect(Rw, Tw, T, x, μ/ρ, k*Pr^cal.n,cal.n))
end

struct GlassbeadCorrect{U<:Real,Fluid} <: AbstractAnemCorrect
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

    
function GlassbeadCorrect(c::GlassbeadCorrect;
                          T=reftemp(c), Rw=resistance(c),
                          Tw=temperature(c), x=fluid(c), P=101325.0)
    ρ = density(x, T, P)
    μ = viscosity(x, T, P)
    cₚ = specheat(x, T, P)
    k = heatcond(x, T, P)
    
    Pr = cₚ * μ / k
    ν  = μ / ρ
    
    meas = GlassbeadCorrect(Rw, Tw, T, x, ν,
                            k*Pr^cal.n, cal.n, cal.β, cal.c1, cal.c2)
    
    GlassbeadCorrect(Rw, Tw, T, c.fluid, c.ν, c.ϕ, c.n,
                     c.β, c.c1, c.c2)
end


function mf58correct(Rsens ;R=1e2, T=298.15, P=101_325.0, n=1/3,
                     N=2, beta=150.0, fluid=Air)
    Tw = temperature(Rsens, R)
    Rw = R
    ρ = density(fluid, T, P)
    μ = viscosity(fluid, T, P)
    kₐ = heatcond(fluid, T, P)
    cₚ = specheat(fluid, T, P)
    
    ν = μ / ρ
    Pr = cₚ * μ / k
    ϕ = kₐ * Pr^n

    D = 2.0e-3 # Diameter of thew glass body
    L = 4.0e-3 # Length of the glass body
    A = π*D*L  # Approximate surface area of the glass body
    c1 = A / D
    
    d = 0.5e-3 # Diameter of the electric steel leads
    P = π*d  # Perimiter of the leads
    Af = π*d^2/4 # Cross section area of the leads
    kf = 30.0 #W/mK # Heat conductivity of steel
    
    γ = sqrt(d/D) # Convection coefficient factor
    
    β = beta  # Ts = Tw - β⋅Q̇
    c2 = N * sqrt(γ*kf*Af*P/D)
    return GlassbeadCorrect(Rw, Tw, T, fluid, ν, ϕ, n, β, c1, c2)
end   

function anemcorrect(E,cal::GlassbeadCorrect, T, P, x, Rw, Tw)
    
    ρ = density(x, T, P)
    μ = viscosity(x, T, P)
    cₚ = specheat(x, T, P)
    k = heatcond(x, T, P)
    
    Pr = cₚ * μ / k
    ν  = μ / ρ
    
    meas = GlassbeadCorrect(Rw, Tw, T, x, ν,
                            k*Pr^cal.n, cal.n, cal.β, cal.c1, cal.c2)
    return anemcorrect(E, cal, meas)
end

function anemcorrect(E, cal::GlassbeadCorrect, meas::GlassbeadCorrect)
    
    Rwc = resistance(cal)
    Twc = temperature(cal)
    Tac = reftemp(cal)

    Rw = resistance(meas)
    Tw = temperature(meas)
    Ta = reftemp(meas)

    β = cal.β
    ϕc = cal.ϕ
    ϕ = meas.ϕ
    
    Y = E*E / (Rw * (Tw - Ta))
    X = Y / (1-β*Y)

    c1 = cal.c1
    c2 = cal.c2

    u = (-c2 + sqrt(c2*c2 + 4*c1*X)) / (2c1)

    fc = u*u / ϕ
    
    Xc = c1*fc*ϕc + c2*sqrt(fc*ϕc)
    Yc = Xc / (1 + β*Xc)
    
    return E*sqrt(Yc * Rwc * (Twc - Tac))
    
end

