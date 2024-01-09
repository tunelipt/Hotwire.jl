# Temperature/Pressure/etc correction


"""
    `AbstractAnemCorrect`

Abstract class that corrects the output of thermal anemometers due to changes
in fluid properties and electronic configurations

The output of a thermal anemometer is directly related to the heat transfer from a 
"""
abstract type AbstractAnemCorrect end



struct CorrFactor{U,V,W}
    "Correction factor"
    f::U
    "Kinematic viscosity"
    nu::V
    "Output Voltage"
    E::W
end
Base.broadcastable(f::CorrFactor) = Ref(f)

kinvisc(f::CorrFactor) = f.nu
corrfactor(f::CorrFactor) = f.f
voltage(f::CorrFactor) = f.E

function (fc::CorrFactor)(Eo, gain, offset)
    Ew = Eo/gain + offset
    return (Ew*fc.f - offset)*gain
end

function (fc::CorrFactor)(Eo, gain, offset, idx::Integer)
    Ew = Eo/gain + offset
    return (Ew*fc.f[i] - offset[i])*gain[i]
end


struct TempCorrect{U,Fluid} <: AbstractAnemCorrect
    "Fluid temperature"
    T::U
    "Operating resistance"
    Rw::U
    "Operating resistance temperature"
    Tw::U
    "Operating pressure"
    P::U
    "Operating fluid"
    fluid::Fluid
    "Kinematic viscosity of the fluid"
    ν::U
end

Base.broadcastable(cal::TempCorrect) = Ref(cal)


resistance(c::AbstractAnemCorrect) = c.Rw
temperature(c::AbstractAnemCorrect) = c.Tw
reftemp(c::AbstractAnemCorrect) = c.T
kinvisc(c::AbstractAnemCorrect) = c.ν
pressure(c::AbstractAnemCorrect) = c.P
fluid(c::AbstractAnemCorrect) = c.fluid

"""
`TempCorrect(T, Rw, Tw, P, fluid, ν)`
`TempCorrect(c::TempCorrect; T=reftemp(c), P=Rw=resistance(c), Tw=temperature(c))`

Creates a `TempCorrect` object. The `TempCorrect` struct stores information on
operating conditions of a thermal anemometer. The copy constructor creates a new
`TempCorrect` object changing individual fields.
"""
function TempCorrect(c::TempCorrect, T, P, fluid, Rw, Tw)
    
    ν = kinvisc(fluid, (T+Tw)/2, P)
    TempCorrect(T,Rw,Tw,P,fluid,ν)
end

TempCorrect(T, P, fluid, Rw, Tw) = TempCorrect(T, Rw, Tw, P, fluid,
                                               kinvisc(fluid, (T+Tw)/2, P))

"""
`tempcorrect(c::AC, op::AC) where {AC<:AbstractAnemCorrect}`
`tempcorrect(c::AbstractAnemCorrect, T, Rw, Tw)`
`tempcorrect(Tc, Rwc, Twc, T, Rw, Tw)`

Calculates the correction factor due to changes in fluid and sensor temperature
(and resistance).
Given a thermal anemometer with operating conditions specified by `c`,
 and the output of the thermal anemometer at different operating conditions given by
`op::TempCorrect` or by fluid  temperature `T`, operating resistance `Rw` and operating temperature `Tw`.
"""
tempcorrect(Tc,Rwc,Twc, T,Rw,Tw) = sqrt( Rwc/Rw  * (Twc-Tc)/(Tw-T) )

tempcorrect(cal::AbstractAnemCorrect, op::AbstractAnemCorrect) =
    sqrt(resistance(cal)/resistance(op) *
    (temperature(cal) - reftemp(cal)) / (temperature(op) - reftemp(op)))

tempcorrect(cal::AbstractAnemCorrect, T, Rw, Tw) =
    tempcorrect(reftemp(cal), resistance(cal), temperature(cal),
                T, Rw, Tw)

correctmodel(mc::TempCorrect, T, P, fluid, Rw, Tw) =
    TempCorrect(mc, T, P, fluid, Rw, Tw)


"""
`correct(E, mc_cal, mc)`
`correct(E, mc_cal, T, P, fluid, Rw, Tw)`
Corrects the output of a thermal anemometer with respect to varying operating
conditions such as electronics, fluid composition, temperature and pressure.

## Arguments

 * `E`: Voltage accross the sensor resistance
 * `mc_cal`: Correction model [`AbstractAnemCorrect`](@ref) at calibration conditions
 * `mc`: Correction model [`AbstractAnemCorrect`](@ref) at operating conditions
 * `T`: Fluid temperature at operating conditions
 * `P`: Fluid pressure at operating conditions
 * `fluid`: Fluid
 * `Rw`: Sensor resistance at operating conditions
 * `Tw`: Sensor temperature at operating conditions

## Returns

Corrected anemometer output at calibration conditions.

## See also

[`TempCorrect`](@ref), [`WireCorrect`](@ref), [`GlassbeadCorrect`](@ref)



"""
function correct(E, cal::TempCorrect, T, P, fluid, Rw, Tw)
    f = tempcorrect(cal, T, Rw, Tw)
    ν = kinvisc(fluid, (T+Tw)/2, P)
    return CorrFactor(f, ν, f*E)
end

function correct(E, cal::TempCorrect, op::TempCorrect)
    f = tempcorrect(cal, op) 
    ν = kinvisc(op)
    return CorrFactor(f, ν, f*E)
end


struct WireCorrect{U<:Real,Fluid} <: AbstractAnemCorrect
    "Fluid temperature"
    T::U
    "Operating resistance"
    Rw::U
    "Operating resistance temperature"
    Tw::U
    "Operating pressure"
    P::U
    "Operating fluid"
    fluid::Fluid
    "Kinematic viscosity of the fluid"
    ν::U
    "Calibration k⋅Prⁿ"
    ϕ::U
    "Prandtl number exponent"
    n::U
end



"""
`WireCorrect(T, P, fluid, Rw, Tw, ν, ϕ, n)`
`WireCorrect(c::WireCorrect, T, P, fluid, Rw, Tw)`
`WireCorrect(T, P, fluid, Rw, Tw, n=1/3)`

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
function WireCorrect(T, P, fluid, Rw, Tw, n=1/3)
    #function WireCorrect(c::WireCorrect, T, P, fluid, Rw, Tw)
    Tm = (T+Tw)/2  # Film temperature
    ρ = density(fluid, Tm, P)
    μ = viscosity(fluid, Tm, P)
    cₚ = specheat(fluid, Tm, P)
    k = heatcond(fluid, Tm, P)

    Pr = cₚ * μ / k
    ν  = μ / ρ
    ϕ = k * Pr^n
    
    WireCorrect(T, Rw, Tw, P, fluid, ν, ϕ, n)
end

WireCorrect(c::WireCorrect, T, P, fluid, Rw, Tw) =
    WireCorrect(T, P, fluid, Rw, Tw, c.n)

correctmodel(mc::WireCorrect, T, P, fluid, Rw, Tw) =
    WireCorrect(T, P, fluid, Rw, Tw, mc.n) 

function correct(E, mc_cal::WireCorrect, mc::WireCorrect) 
    
    f = tempcorrect(mc_cal, mc) * sqrt(mc_cal.ϕ/mc.ϕ) 
    
    return CorrFactor(f, kinvisc(mc), E*f)
end

function correct(E, mc_cal::WireCorrect, T, P, fluid, Rw, Tw)
    mc = WireCorrect(mc_cal, T, P, fluid, Rw, Tw)
    return correct(E, mc_cal, mc)
end


struct GlassbeadCorrect{U<:Real,Fluid} <: AbstractAnemCorrect
    "Fluid temperature"
    T::U
    "Operating resistance"
    Rw::U
    "Operating resistance temperature"
    Tw::U
    "Operating pressure"
    P::U
    "Operating fluid"
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

"""
`GlassbeadCorrect(T, Rw, Tw, P, fluid, ν, ϕ, n, β, c1, c2)`
`GlassbeadCorrect(c::GlassbeadCorrect, T, P, fluid, Rw, Tw)`

Implements a correction model where the resistor element is surrounded by
by an insulating shell and the heat generated by the resistance element
is transfered to the fluid by

 * The outer surface of the insulating shell
 * The wire leads working as fins

For details in the correction model, see the manual. As an actual example
using `GlassbeadCorrect`, see [`mf58correct`](@ref) for 

## Arguments
 * `ν` Kinematic viscosity of the fluid
 * `ϕ` ``k⋅Prⁿ`` of the fluid
 * `n` Exponent of ``Pr``number when computing `ϕ`
 * `β` Thermal resistivity of the insulating shell such that ``T_s = T_w - β⋅Q̇_w``
 * `c1` and `c2` coefficients such that ``Q̇_w = \\left[c1⋅f(Re)ϕ + c2⋅z\\sqrt{f(Re)ϕ}\\right]⋅(T_s - T_a)``  (notice that it is ``T_s`` and not ``T_w``!)

When only the fluid properties change, a new `GlassbeadCorrect` object can be created
maintaining everything else and the following arguments are used
 * `c` A `GlassbeadCorrct` object that will be modified. Usually the calibration value.
 * `tc` [`TempCorrect`](@ref) object specifying new operating conditions
 * `fluid` Object used to calculate new relevant properties of the fluid
 * `P`  Atmospheric pressure in ``Pa``

"""
function GlassbeadCorrect(c::GlassbeadCorrect, T, P, fluid, Rw, Tw)
    
    Tm = (T+Tw) / 2
    ρ = density(fluid, Tm, P)
    μ = viscosity(fluid, Tm, P)
    cₚ = specheat(fluid, Tm, P)
    k = heatcond(fluid, Tm, P)
    
    Pr = cₚ * μ / k
    ν  = μ / ρ
    
    return GlassbeadCorrect(T, Rw, Tw, P, fluid, ν, k*Pr^c.n, c.n, c.β, c.c1, c.c2)
    
end

correctmodel(mc::GlassbeadCorrect, T, P, fluid, Rw, Tw) =
    GlassbeadCorrect(mc, T, P, fluid, Rw, Tw) 

"""
`mf58correct(T, P, fluid, Rw, Tw; n=1/3, q=0.4,
             N=2, beta=150.0, fluid=Air,
             D=2e-3, L=4e-3, d=0.5e-3, kf=30.0)`

Creates a [`GlassbeadCorrect`](@ref) for thermistors of type
MF-58.

These thermistors can have a wide variation in nominal resistance but they have
two leads made of electric steel with a diameter of 0.5 mm. The body itself is made
of glass and has a diameter of 2 mm and length of 4 mm (tapered at the ends).

From finite element simulations and guessing the exact thermal conductivity of the glass
the value for the thermal resistance ``β`` is around 150 ``K/W``. Further studies
 (numerical and preferentially experimental) should be carried out to estimate more
accurate values for ``β``

Of course, the MF-58 thermistor have two wire leads but depending on probe construction,
both leads might be insulated so the heat they loose to the fluid is negligible.

The convective heat transfer coefficient in the glass body is different from the value
observed for the leads. In this model, we assume that ``h_f`` (leads) is a fixed
fraction (``γ``) of the value for the body ``h_a``. To estimate this fraction, we assume
that

``Nu(Re,Pr) = a⋅Re^q ⋅ Prⁿ``

With this correlation, we have

``γ = h_f / h_a = (d / D)^q ``

## Arguments
 * `T` Fluid temperature
 * `Rw` Resistance of the resistance sensor
 * `Tw` Temperature of the resistance sensor
 * `P` Atmospheric pressure in ``Pa``
 * `fluid` Fluid 
 * `n` Assumed exponent for ``Nu(Re,Pr) = f(Re)⋅Prⁿ`` correlation
 * `q` Re exponent in assumed correlation for ``Nu`` see above explanation
 * `N` Number of _uninsulated_ leads leaving the body
 * `beta` Value the thermal resistance of the glass insulator in ``K/W``
 * `fluid` Fluid used for calculating the properties
 * `D` Diameter in ``mm`` of the glass bead 
 * `L` Length in ``mm`` of the glass bead
 * `d` Diameter in ``mm`` of the leads
 * `kf` Thermal conductivity of the leads

## Return value

A [`GlassbeadCorrect`](@ref) object
"""
function mf58correct(T, P, fluid, Rw, Tw; n=1/3, q=0.4,
                     N=2, beta=150.0, D=2e-3, L=4e-3, d=0.5e-3, kf=30.0)
    Tm = (Tw + T) / 2
    ρ = density(fluid, Tm, P)
    μ = viscosity(fluid, Tm, P)
    kₐ = heatcond(fluid, Tm, P)
    cₚ = specheat(fluid, Tm, P)
    
    ν = μ / ρ
    Pr = cₚ * μ / kₐ
    ϕ = kₐ * Pr^n

    #D = 2.0e-3 # Diameter of thew glass body
    #L = 4.0e-3 # Length of the glass body
    A = π*D*L  # Approximate surface area of the glass body
    c1 = A / D
    
    #d = 0.5e-3 # Diameter of the electric steel leads
    Pf = π*d  # Perimiter of the leads
    Af = π*d^2/4 # Cross section area of the leads
    #kf = 30.0 #W/mK # Heat conductivity of steel
    
    γ = sqrt(d/D) # Convection coefficient factor
    
    β = one(ϕ)*beta  # Ts = Tw - β⋅Q̇
    c2 = N * sqrt(γ*kf*Af*Pf/D)
    return GlassbeadCorrect(T, Rw, Tw, P, fluid, ν, ϕ, n, β, c1, c2)
end   


function correct(E, mc_cal::GlassbeadCorrect, mc::GlassbeadCorrect)
    
    Rwc = resistance(mc_cal)
    Twc = temperature(mc_cal)
    Tac = reftemp(mc_cal)

    Rw = resistance(mc)
    Tw = temperature(mc)
    Ta = reftemp(mc)

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
    E1 = sqrt(Yc * Rwc * (Twc - Tac))
    return CorrFactor(E1/E, kinvisc(mc), E1)
    
end

correct(E, mc_cal::GlassbeadCorrect, T, P, fluid, Rw, Tw) = 
    correct(E, mc_cal, correctmodel(mc_cal, T, P, fluid, Rw, Tw))

