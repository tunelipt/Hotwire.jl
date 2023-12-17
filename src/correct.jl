# Temperature/Pressure/etc correction





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
anemcorrectfactor(cal::TempCorrect, op::TempCorrect) =
    sqrt(resistance(cal)/resistance(op) *
    (temperature(cal) - reftemp(cal)) / (temperature(op) - reftemp(op)))

anemcorrectfactor(cal::TempCorrect; T=reftemp(cal), Rw=resistance(cal),
                  Tw=temperature(cal)) = anemcorrectfactor(cal, TempCorrect(T,Rw,Tw))

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

kinvisc(c::AbstractAnemCorrect) = c.ν

"""
`anemcorrect(E, tc_cal, mc_cal, tc, mc)`

Corrects the output of a thermal anemometer with respect to varying operating
conditions such as electronics, fluid composition, temperature and pressure.

## Arguments

 * `E`: Voltage accross the sensor resistance
 * `tc_cal`: [`TempCorrect`](@ref) at calibration conditions
 * `mc_cal`: A correction model for operating conditions at calibration conditions
 * `tc`: [`TempCorrect`](@ref) at operating conditions
 * `mc`: A correction model for operating conditions at operating conditions

## Returns

Corrected anemometer output at calibration conditions.

## See also

[`TempCorrect`](@ref), [`WireCorrect`](@ref), [`GlassbeadCorrect`](@ref)



"""
anemcorrect(E, tc_cal::TempCorrect, mc_cal, tc::TempCorrect, mc) =
    E*anemcorrectfactor(tc_cal, tc)

"""
`TOnlyCorrect(ν)`

No correction model, only [`TempCorrect`](@ref) will be used to correct anemometer output.
Only a change in kinematic viscosity will lead to a change in flow velocity.
"""
struct TOnlyCorrect{U} <: AbstractAnemCorrect
    "Kinematic viscosity"
    ν::U
end

function TOnlyCorrect(c::TOnlyCorrect, tc::TempCorrect, fluid, P=101325.0)
    T = reftemp(tc)
    ν = kinvisc(fluid, T, P)
    return TOnlyCorrect(ν)
end

function anemcorrect(E, tc_cal, mc_cal::AC, T, Tw, Rw, fluid, Pa) where{AC}
    tc = TempCorrect(T,Rw,Tw)
    mc = AC(mc_cal, tc, fluid, Pa)
    return anemcorrect(E, tc_cal, mc_cal, tc, mc)
end
#WireCorrect(c::WireCorrect, tc::TempCorrect, fluid, P=101325.0) =

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
`WireCorrect(tc::TempCorrect, fluid, P, n)`
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
function WireCorrect(tc::TempCorrect, fluid, P=101325.0, n=1/3)
    T = reftemp(tc)
    Tw = temperature(tc)
    Tm = (T+Tw)/2  # Film temperature
    ρ = density(fluid, Tm, P)
    μ = viscosity(fluid, Tm, P)
    cₚ = specheat(fluid, Tm, P)
    k = heatcond(fluid, Tm, P)

    Pr = cₚ * μ / k
    ν  = μ / ρ
    ϕ = k * Pr^n
    
    WireCorrect(ν, ϕ, n)
end
WireCorrect(c::WireCorrect, tc::TempCorrect, fluid, P=101325.0) =
    WireCorrect(tc, fluid, P, c.n)


function anemcorrect(E, tc_cal::TempCorrect, mc_cal::WireCorrect,
                     tc::TempCorrect, mc::WireCorrect) 
    
    f = anemcorrectfactor(tc_cal, tc)
    return E*sqrt(mc_cal.ϕ/mc.ϕ) * f
end

struct GlassbeadCorrect{U<:Real} <: AbstractAnemCorrect
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
`GlassbeadCorrect(c::GlassbeadCorrect, tc::TempCorrect, fluid, P)`

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
function GlassbeadCorrect(c::GlassbeadCorrect{U}, tc::TempCorrect,
                          fluid, P=101325.0) where {U}

    T = reftemp(tc)
    Tw = temperature(tc)
    Tm = (T+Tw) / 2
    ρ = density(fluid, Tm, P)
    μ = viscosity(fluid, Tm, P)
    cₚ = specheat(fluid, Tm, P)
    k = heatcond(fluid, Tm, P)
    
    Pr = cₚ * μ / k
    ν  = μ / ρ
    
    return GlassbeadCorrect{U}(U(ν), U(k*Pr^c.n), U(c.n), U(c.β), U(c.c1), U(c.c2))
    
end


"""
`mf58correct(tc::TempCorrect; P=101_325.0, n=1/3, q=0.4,
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
 * `tc`  [`TempCorrect`](@ref) object specifying the operating conditions of the anemometer
 * `P` Atmospheric pressure in ``Pa``
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
function mf58correct(tc::TempCorrect, fluid=AIR, P=101_325.0; n=1/3, q=0.4,
                     N=2, beta=150.0, D=2e-3, L=4e-3, d=0.5e-3, kf=30.0)
    Tw = temperature(tc)
    Rw = resistance(tc)
    T = reftemp(tc)
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
    return GlassbeadCorrect(ν, ϕ, n, β, c1, c2)
end   


function anemcorrect(E, tc_cal::TempCorrect, mc_cal::GlassbeadCorrect,
                     tc::TempCorrect, mc::GlassbeadCorrect) 

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
    
    return sqrt(Yc * Rwc * (Twc - Tac))
    
end

