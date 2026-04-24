"""
Abstract type for handling calibrations and corrections

The calibration of a thermal anemometer is basically a table of output voltages (E)
for different velocities (U). This is carried out at a specific temperature, pressure
and fluid composition.

If during use of the anemometer, everything is the same, this calibration table can be
used directly. But, what happens when something changes? The most important change
is when fluid temperature changes. The heat transfer between the surface of the sensor
and the fluid changes and so does the output of the anemometer that should be corrected
to the calibration conditions.

Typical changes are temperature, pressure, fluid and operating conditions of the sensor
(electronics).

Another important aspect is the operation of the sensor. The most common methods of
operating a thermal anemometer is

 * Constant temperature, where the temperature of the resistive element is kept constant
 * Constant current, where the current flowing through the resistive element is kept
constant

An `AbstractAnemCalibr` handles the calibration and correction. For now, constant
temperature (CTA) and constant current (CCA) are handled.

CTA is the simpler case since the temperature of the resistive element is constant.
This temperature corresponds to a constant resistance. But in the case of CCA, since
the current is constant, the temperature (and therefore resistance) will depend on
the heat transfer and fluid velocity.

To use the same approach for both CTA and CCA anemometers, the calibration conditions
are corrected to a standard calibration resistance (and calibration temperature).

The interface for `AbstractAnemCalibr` objects is:

```julia

CalibrationMethod(R::RT,
                  E::AbstractVector,
                  U::AbstractVector, Rw,
                  T, P, makefitfun;
                  fluid=AIR, params...) where {RT<:AbstractResistor}
```

The input arguments are
 * `R` an [`AbstractResistor`](@ref) object ([`Resistor`](@ref) or [`Thermistor`](@ref))
 * `E` `AbstractVector` containing the voltage output of the sensor during calibration
 * `U` `AbstractVector` containing the calibration velocity
 * `Rw` Op. calibration resistance. Constant value for CTA but a vector otherwise.
 * `T` Calibration temperature - could be a scalar or a vector
 * `P` Calibration pressure - could be a scalar or a vector
 * `makefitfun` a function the creates a curve fit for calibration data.
 
Key word arguments
 * `fluid` Object that allows the calculation of thermodynamic and transport properties
 * `params...` model specific parameters.

#### Basic approach
The heat transfer from a heated element is given by Newton's law of cooling:

```math
Q̇ = h(U, ...) (Tw - Ta)
```

But ``Q̇`` is dependent on the electric properties of the self heated resistor:
```math
Q̇ = E⋅I = E² / Rw
```

The convection coefficient ``h``, from dimensional analysis is dependent on the Reynolds number, therefor
```math
E² / (Rw⋅(Tw-Ta)⋅stuff...) = f(U/ν)
```
where stuff is what the different models will modify. In the simplest case where fluid
properties are constant, stuff is just 1.

inverting this relationship we have
```math
U/ν = g( E² / (Rw⋅(Tw-Ta)⋅stuff...) )
```

#### Curve fit

The curve fitting scheme should fit the above equations. Different possibilities exist but it helps to rememember that in general

```math
Q̇ ∝ √U
```

So something like
```math
(U/ν)ᵅ = a₀ + a₁E² + a₂E⁴ + a₃E⁶ + …
```

works really well.
The appropriate curve fitting procedure depends on the exact model and can not be
specified here. At first thought it should be a curve (possibly a polynomial)
of type ``U = U(E)`` but different models use different approaches. [`TempCalibr`](@ref)
uses this but [`HWCalibr`](@ref) fits different variables.

The fluid is specified by an object that needs to provide the following methods:

 * [`heatcond`](@ref) Heat conductivity
 * [`viscosity`](@ref) Dynamic viscosity
 * [`density`](@ref) Fluid density
 * [`specheat`](@ref) Specific heat
 * [`prandtl`](@ref) Prandtl number
 * [`kinvisc`](@ref) Kinematic viscosity.

Any of the above properties should be used as `property(fluid, T, P)`.

For constant properties fluids, there is the [`ConstPropFluid`](@ref) type.
The [`IdealGas`](@ref) where the properties of the fluid are calculated for Ideal
gases from the molecular weight, `Cp(T)`, `μ(T)`, `k(T)`. A few common gases are
implemented such as `AIR`, `NITROGEN`, `OXYGEN`, `HELIUM`, `C3H8`, `CO2` and
`HYDROGEN`.

"""
abstract type AbstractAnemCalibr end

"Calibration pressure"
pressure(cal::AbstractAnemCalibr) = cal.P
calpress(cal::AbstractAnemCalibr) = cal.P

"Calibration temperature"
temperature(cal::AbstractAnemCalibr) = cal.T
caltemp(cal::AbstractAnemCalibr) = cal.T

"Calibration fluid"
fluid(cal::AbstractAnemCalibr) = cal.fluid
calfluid(cal::AbstractAnemCalibr) = cal.fluid

"""
`makecaltable(E, U, T, P)`

Returns a matrix with calibration conditions.

The temperature and pressure can be scalars or vectors.
"""
function makecaltable(E::AbstractVector{X}, U::AbstractVector{Y}, T, P) where {X,Y}
    N = length(E)
    
    @assert size(E,1) == size(U,1)

    caltab = zeros(promote_type(X,Y), N, 4)

    caltab[:,1] .= E
    caltab[:,2] .= U
    caltab[:,3] .= T
    caltab[:,4] .= P

    Tm = mean(T)
    Pm = mean(P)
    
    return Tm, Pm, caltab
end


struct TempCalibr{X,Fluid,Fit} <: AbstractAnemCalibr
    "Calibration temperature"
    T::X
    "Calibration pressure"
    P::X
    "Calibration fluid"
    fluid::Fluid
    "Calibration curve fit"
    fit::Fit
end
Base.broadcastable(calibr::TempCalibr) = Ref(calibr)

"""
`TempCalibr(R, Ec, Uc, Tc, Pc, Rwc, makefitfun;
            Rw=nothing, T=nothing, P=nothing, fluid=AIR)`

See [`AbstratcAnemCalibr`](@ref) for a general description of the interface.

The `TempCalibr` model assumes that a same fluid with constant fluid properties is
used.

The output voltage of the anemometer is given by
```math
E² = Rw ⋅ h(U) ⋅ (Tw - Ta)
```

If properties are constant, then the output voltage in calibration conditions is
given by:

```math
Ec² = E² ⋅ Rwc/Rw  ⋅ (Twc - Tac)/(Tw - Ta)
```

#### Curve fit

This model fits directly ``U=U(E)``. A simple, not the most accurate but robust curve
fit is the generalized King's law. It fits

```math
E² = A + B⋅Uⁿ
```

This is implemented in [`KingLaw`](@ref).

Simple n-degree polynomials are also commonly used but they tend to oscillate. A better and more accurate possibility is what we are calling King's Polynomial:

```math
Uⁿ = a₀ + a₁ E² + a₂ E⁴ + ... + aₖ (E²)ᵏ
```
The `n` parameter can be fitted or imposed. A value of 0.45 usually gets good results.

This fit is implemented in [`KingPoly`](@ref).

"""
function TempCalibr(R::RT,
                    U::AbstractVector,
                    E::AbstractVector, Rw,
                    T, P, makefitfun;
                    fluid=AIR) where {RT<:AbstractResistor}
    
    @assert size(U) == size(E) 
    Tw = temperature.(R, Rw)
    
    ΔT = Tw .- T
    
    fit = makefitfun(E.^2 ./ (ΔT .* Rw), U)
    
    
    return TempCalibr(mean(T), mean(P), fluid, fit)
    
end

function hwcorrect(cal::TempCalibr, R::RT, Rw, T, P,
                   fluid) where {RT<:AbstractResistor}
    Tw = temperature(R, Rw)
    ΔT = Tw - T
    return Rw*ΔT, 1
end

"""
`velocity(cal, E, Rw, T, P, fluid, Rw)`
`velocity(cal, E, Rw; T, P, fluid, Rw)`

Calculate the velocity given anemometer output `E` and flow conditions `T`, `P`, with fluid `fluid` and operating at resistance `Rw`.

"""    
function velocity(cal::TempCalibr, R::RT, E, Rw, T, P,
                  fluid) where {RT<:AbstractResistor}
    RΔT, ν = hwcorrect(cal, R, Rw, T, P, fluid)
    return ν * cal.fit(E*E/RΔT) 
end



struct HWCalibr{X<:AbstractFloat,Fluid,Fit} <: AbstractAnemCalibr
    "Calibration temperature"
    T::X
    "Calibration pressure"
    P::X
    "Calibration fluid"
    fluid::Fluid
    "Calibration curve fit"
    fit::Fit
    "Prandtl number exponent in Nusselt"
    n::X
    "Mean surface temperature factor"
    theta::X
end
Base.broadcastable(calibr::HWCalibr) = Ref(calibr)

"""
`HWCalibr(R, Ec, Uc, Tc, Pc, Rwc, makefitfun;
            Rw=nothing, T=nothing, P=nothing, fluid=AIR,
            n=0.3, theta=0.5)`

This model is inspired by Hultmark and Smits (2010) paper.

### Arguments
 * `R` an [`AbstractResistor`](@ref) object ([`Resistor`](@ref) or [`Thermistor`](@ref))
 * `Ec` `AbstractVector` containing the voltage output of the sensor during calibration
 * `Uc` `AbstractVector` containing the calibration velocity
 * `Tc` Calibration temperature - could be a scalar or a vector
 * `Pc` Calibration pressure - could be a scalar or a vector
 * `Rwc` Op. calibration resistance. Constant value for CTA but a vector otherwise.
 * `makefitfun` a function the creates a curve fit for calibration data.

Key word arguments
 * `Rw` specific calibration operating resistance. If not given use mean of `Rwc`
 * `T` specific calibration temperature. If not given, use mean of `Tc`
 * `P` specific calibration pressure. If not given, use mean of `Pc`
 * `fluid` Object that allows the calculation of thermodynamic and transport properties
 * `n` Pr exponent for calculating Nusserlt
 * `theta` Fraction of (Tw - Ta) used for calculating transport and thermodynamic properties of the fluid.


### Model 
The heat transfer from the sensor surface is given by
```math
Q̇ = E²/Rw = h(U,...)⋅A⋅(Ts - Ta)
```
The convection coefficient ``h`` is calculated from the Nusselt number:
```math
Nu = hL/k = Nu(Re,Pr)
```
where ``k`` is the thermal conductivity of the fluid, Re is the Reynolds number and Pr is the Prandtl number. This model assumes that the Nu can be calculated as
```math
Nu(Re,Pr) = f₁(Re) Prⁿ
```
where n≈0.3 is the first model coefficient that should be provided. With such an approximation,
```math
hA = k/L⋅f₁(Re) Prⁿ A = k⋅Prⁿ⋅f₂(U/ν)
```

Given the calibration table, this model uses the following equation for fitting the data:
```math
E²/[k Prⁿ Rw (Tw - Ta)] = f₂(U/ν)
```
or,
```math
U/ν = f₃( E²/[k Prⁿ Rw (Tw - Ta)] )
```

What temperature should be used to calculata transport and thermodynamic properties of the fluid? Traditionally, the film temperature is used, the mean between surface temperature and fluid temperature. But in case the surface temperature is not equal to Tw, a possible situation if cheap thermistors are used, then some other temperature should be used. The second parameter is θ, so that the surface temperature is approximately equal to:

```math
Ts = Ta + θ(Tw-Ta)
```

When ``Ts=Tw``, θ=1/2. This value should be independently estimated, either through numerical simulations or direct measurement. Keep in mind that θ is only used to calculate the fluid properties.

### Curve Fit

The model will fit

```math
U/ν = f₃( E²/[k Prⁿ Rw (Tw - Ta)] )
```
An often used correlation for Nusselt is
```math
Nu(Re,Pr) = a Reᵐ Prⁿ
```
so it is reasonable to approximate the curve above as

```math
(U/ν)ᵐ = a₀ + a₁E + a₂E² + a₃E³ + ... + aₖEᵏ
```

"""
function HWCalibr(R::RT,
                  U::AbstractVector,
                  E::AbstractVector, Rw,
                  T, P, makefitfun;
                  fluid=AIR, n=0.3, theta=0.5) where {RT<:AbstractResistor}
    
    Tw = temperature.(R, Rw)
    
    ΔT = Tw .- T
    
    # Film temperature
    Tf  = T   .+ theta  .* ΔT  # Calibration reference conditions
    phi   = heatcond.(fluid, Tf,  P) .* prandtl.(fluid, Tf, P) .^ n

    # Reynolds number
    nu = kinvisc.(fluid, Tf, Pc)
    Re = U ./ nu

    # E²/ϕΔT
    Ex = E .* E ./ (phi .* Rw .* ΔT)

    fit = makefitfun(Ex, Rec)
    return HWCalibr(mean(T), mean(P), fluid, fit, n, theta)
    
end

function hwcorrect(cal::HWCalibr, R::RT, Rw, T, P,
                   fluid) where {RT<:AbstractResistor}
    Tw = temperature(R, Rw)
    ΔT = Tw - T
    Tf = T + θ*ΔT
    
    k = heatcond(fluid, Tf, P)
    Pr = prandtl(fluid, Tf, P)
    phi = k * Pr^cal.n
    ν = kinvisc(fluid, Tf, P)
    return Rw*phi*ΔT, ν
end


function velocity(cal::HWCalibr, R::RT, E, Rw, T, P,
                  fluid) where {RT<:AbstractResistor}
    RϕΔT, ν = hwcorrect(cal, R, Rw, T, P, fluid)
    return ν * cal.fit(E*E/RϕΔT)    
end



