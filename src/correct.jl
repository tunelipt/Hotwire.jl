# Temperature/Pressure/etc correction

export AbstractAnemCorrect, TempCorrect, WireCorrect, GlassbeadCorrect, TOnlyCorrect
export mf58correct, anemcorrect, anemcorrectfactor
    




struct TempCorrect{U<:Real,V<:Real,W<:Real}
    "Fluid temperature"
    T::U
    "Operating resistance"
    Rw::V
    "Operating resistance temperature"
    Tw::W
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
`anemcorrectfactor(c::TempCorrect, op::TempCorrect)`
`anemcorrectfactor(c::TempCorrect,T=reftemp(cal),Rw=resistance(cal),Tw=temperature(cal))`

Calculates the correction factor due to changes in operating conditions.
Given a thermal anemometer with operating conditions specified by `cal`,
 and the output of the thermal anemometer at different operating conditions given by
`op::TempCorrect` or by fluid  temperature `T`, operating resistance `Rw` and operating temperature `Tw`.
"""
anemcorrectfactor(cal::TemCorrect, op::TempCorrect) =
    sqrt(resistance(cal)/resistance(op) *
    (temperature(cal) - reftemp(cal)) / (temperature(op) - reftemp(op)))

anemcorrectfactor(cal::TempCorrect; T=reftemp(cal), Rw=resistance(cal),
                  Tw=temperature(cal)) = anemcorrectfactor(cal, TempCorr(T,Rw,Tw))

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

kinvisc(c::AbstractAnemCorrec) = c.ν

"""
`anemcorrect(E, tc_cal, mc_cal, tc, mc)`

Corrects the output of a thermal anemometer with respect to varying operating
conditions such as electronics, fluid composition, temperature and pressure.

## Arguments

 * `E`: Voltage accross the sensor resistanve
 * `tc_cal`: [`TempCorr`](@ref) at calibration conditions
 * `mc_cal`: A correction model for operating conditions at calibration conditions
 * `tc`: [`TempCorr`](@ref) at operating conditions
 * `mc`: A correction model for operating conditions at operating conditions

## Returns

Corrected anemometer output at calibration conditions.

## See also

[`TempCorr`](@ref), [`WireCorrect`](@ref), [`GlassbeadCorrect`](@ref)



"""
anemcorrect(E, tc_cal::TempCorrect, mc_cal, tc::TempCorrect, mc) =
    E*anemcorrectfactor(tc_cal, tc)

"""
`TOnlyCorrect(ν)`

No correction model, only [`TempCorr`](@ref) will be used to correct anemometer output.
Only a change in kinematic viscosity will lead to a change in flow velocity.
"""
struct TOnlyCorrec{U} <: AbstractAnemCorrect
    "Kinematic viscosity"
    ν::U
end

function TOnlyCorrect(c::WireCorrect, tc::TempCorrect, fluid, P=101325.0)
    T = reftemp(tc)
    ν = kinvisc(fluid, T, P)
    return TOnlyCorrect(ν)
end


function anemcorrect(E, tc_cal::TempCorrect, mc_cal::TOnlyCorrect,
                     tc::TempCorrect, mc::TOnlyCorrect) 
    
    f = anemcorrectfactor(tc_cal, tc)
    return E * f
end

struct WireCorrect{U<:Real,V<:Real,W<:Real} <: AbstractAnemCorrect
    "Calibration kinematic viscosity"
    ν::U
    "Calibration k⋅Prⁿ"
    ϕ::V
    "Prandtl number exponent"
    n::W
end


"""
`WireCorrect(ν, ϕ, n)`
`WireCorrect(c::WireCorrect, tc::TempCorrect, fluid, P)`

Creates `WireCorrect` object. This object allows the correction for thermal anemometer output when the surface of the resistive element _is_ equal to the temperature of the resitive element.

This correction modelo assumes that the Nusselt number of the heated element is

    ``Nu(Re,Pr) = f(Re)⋅Prⁿ``

## Arguments

 * `ν` Kinematic viscosity of the fluid
 * `ϕ` ``k⋅Prⁿ`` of the fluid where ``k`` is the heat conductivity of the fluid and ``Pr`` its Prandtl number
 * `n` is the exponent of ``Pr`` in the correlation. A value of 1/3 should be good.

In order to conform to the API, the second form build a new `WireCorrect` object from a previous (calibration) `WireCorrect` object under new conditions (fluid, pressure and `TempCorr`)

 * `c` A `WireCorrect` that is built upon
 * `tc` specifies the [`TempCorrect`](@ref)
 * `fluid` specifies the fluid
 * `P` Fluid pressure

"""
function WireCorrect(c::WireCorrect, tc::TempCorrect, fluid, P=101325.0)
    T = reftemp(tc)
    
    ρ = density(fluid, T, P)
    μ = viscosity(fluid, T, P)
    cₚ = specheat(fluid, T, P)
    k = heatcond(fluid, T, P)

    Pr = cₚ * μ / k
    ν  = μ / ρ
    ϕ = k * Pr^c.n
    
    WireCorrect(ν, ϕ, n)
end

function anemcorrect(E, tc_cal::TempCorrect, mc_cal::WireCorrect,
                     tc::TempCorrect, mc::WireCorrect) 
    
    f = anemcorrectfactor(tc_cal, tc)
    return E*sqrt(mc_cal.ϕ/mc.ϕ) * f
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

"""
`GlassbeadCorrect(ν, ϕ, n, β, c1, c2)`
`GlassbeadCorrect(ν, ϕ, n, β, c1, c2)`

"""
function GlassbeadCorrect(c::GlassbeadCorrect, tc::TempCorrect, fluid, P=101325.0)

    T = reftemp(tc)
    ρ = density(x, T, P)
    μ = viscosity(x, T, P)
    cₚ = specheat(x, T, P)
    k = heatcond(x, T, P)
    
    Pr = cₚ * μ / k
    ν  = μ / ρ
    
    return GlassbeadCorrect(ν, k*Pr^c.n, c.n, c.β, c.c1, c.c2)
    
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

function anemcorrect(E, tc_cal::TempCorrect, mc_cal::WireCorrect,
                     tc::TempCorrect, mc::WireCorrect) 

    Rwc = resistance(tc_cal)
    Twc = temperature(tc_cal)
    Tac = reftemp(tc_cal)

    Rw = resistance(tc)
    Tw = temperature(tc)
    Ta = reftemp(tc)

    β = mc_cal.β
    ϕc = mc_cal.ϕ
    ϕ = mc.ϕ
    
    Y = E*E / (Rw * (Tw - Ta))
    X = Y / (1-β*Y)

    c1 = mc_cal.c1
    c2 = mc_cal.c2

    u = (-c2 + sqrt(c2*c2 + 4*c1*X)) / (2c1)

    fc = u*u / ϕ
    
    Xc = c1*fc*ϕc + c2*sqrt(fc*ϕc)
    Yc = Xc / (1 + β*Xc)
    
    return E*sqrt(Yc * Rwc * (Twc - Tac))
    
end

