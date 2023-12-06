# Temperature/Pressure/etc correction

export AbstractAnemCorrect, TempCorrect, WireCorrect, GlassbeadCorrect
export mf58correct
    




struct TempCorrect{U<:Real,V<:Real,W<:Real}
    "Fluid temperature"
    T::W
    "Operating resistance"
    Rw::U
    "Operating resistance temperature"
    Tw::V
end

Base.broadcastable(cal::TempCorrect) = Ref(cal)


resistance(c::TempCorrect) = c.Rw
temperature(c::TempCorrect) = c.Tw
reftemp(c::TempCorrect) = c.T


"""
`TempCorrect(T, Rw, Tw)`
`TempCorrect(c::TempCorrect; T=reftemp(c), Rw=resistance(c), Tw=temperature(c))`

Creates a `TempCorrect` object. The `TempCorrect` struct stores information on
operating conditions of a thermal anemometer. The copy constructor creates a new
`TempCorrect` object changing individual fields.
"""
function TempCorrect(c::TempCorrect;
                     T=reftemp(c), Rw=resistance(c), Tw=temperature(c))
    TempCorrect(Rw, Tw, T)
end

"""
`correctfactor(cal::TempCorrect, op::TempCorrect)`
`correctfactor(cal::TempCorrect, T=reftemp(cal),Rw=resistance(cal),Tw=temperature(cal))`

Calculates the correction factor due to changes in operating conditions.
Given a thermal anemometer with operating conditions specified by `cal`,
 and the output of the thermal anemometer at different operating conditions given by
`op::TempCorrect` or by fluid  temperature `T`, operating resistance `Rw` and operating temperature `Tw`.
"""
anemcorrectfactor(cal::TemCorrect, op::TempCorrect) =
    sqrt(resistance(cal)/resistance(op) *
    (temperature(cal) - reftemp(cal)) / (temperature(op) - reftemp(op)))

anemcorrectfactor(cal::TempCorrect; T=reftemp(cal), Rw=resistance(cal),
                  Tw=temperature(cal)) =
                      sqrt(resistance(cal)/Rw *
                      (temperature(cal)-reftemp(cal)) / (Tw-T))

(cal::TempCorrect)(op::TempCorrect) = anemcorrectfactor(cal, op)
(cal::TempCorrect)(; T=reftemp(cal), Rw=resistance(cal), Tw=temperature(cal)) =
    anemcorrectfactor(cal; T=T, Tw=Tw, Rw=Rw)


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
    "Calibration kinematic viscosity"
    ν::U
    "Calibration k⋅Prⁿ"
    ϕ::U
    "Prandtl number exponent"
    n::U
end

"""
`WireCorrect(c::WireCorrect, cal::TempCorrect, fluid; 
             T=reftemp(cal), Rw=resistance(cal), Tw=temperature(cal),
             P=101325.0)`

Creates an object 
"""
function WireCorrect(c::WireCorrect, cal::TempCorrect, fluid;
                     T=reftemp(cal), Rw=resistance(cal), Tw=temperature(cal),
                     P=101325.0)
    ρ = density(fluid, T, P)
    μ = viscosity(fluid, T, P)
    cₚ = specheat(fluid, T, P)
    k = heatcond(fluid, T, P)

    Pr = cₚ * μ / k
    ν  = μ / ρ
    ϕ = k * Pr^c.n
    
    WireCorrect(Rw, Tw, T, x, ν, ϕ, n)
end

function anemcorrect(E, cmod::WireCorrect, caltemp::TempCorrect,
                     opmod::WireCorrect, optemp::TempCorrect)
    f = anemcorrectfactor(caltemp, optemp)
    return E*sqrt(cmod.ϕ/opmod.ϕ) * f
end

struct GlassbeadCorrect{U<:Real,Fluid} <: AbstractAnemCorrect
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

function anemcorrect(E, cmod::WireCorrect, caltemp::TempCorrect,
                     opmod::WireCorrect, optemp::TempCorrect)
    
    Rwc = resistance(caltemp)
    Twc = temperature(caltemp)
    Tac = reftemp(caltemp)

    Rw = resistance(optemp)
    Tw = temperature(optemp)
    Ta = reftemp(optemp)

    β = cmod.β
    ϕc = cmod.ϕ
    ϕ = opmod.ϕ
    
    Y = E*E / (Rw * (Tw - Ta))
    X = Y / (1-β*Y)

    c1 = cmod.c1
    c2 = cmod.c2

    u = (-c2 + sqrt(c2*c2 + 4*c1*X)) / (2c1)

    fc = u*u / ϕ
    
    Xc = c1*fc*ϕc + c2*sqrt(fc*ϕc)
    Yc = Xc / (1 + β*Xc)
    
    return E*sqrt(Yc * Rwc * (Twc - Tac))
    
end

