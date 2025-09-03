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

CalibrationMethof(R::RT,
                  Ec::AbstractVector,
                  Uc::AbstractVector,
                  Tc, Pc, Rwc, makefitfun;
                  Rw=nothing, T=nothing, P=nothing,
                  fluid=AIR, params...) where {RT<:AbstractResistor}
```

The input arguments are
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
 * `params...` model specific parameters.


#### Curve fit
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

"Calibration operating temperature of the sensor"
temperature(cal::AbstractAnemCalibr) = cal.Tw

"Calibration pressure"
pressure(cal::AbstractAnemCalibr) = cal.P

"Calibration temperature"
caltemp(cal::AbstractAnemCalibr) = cal.T

"Reference temperature"
reftemp(cal::AbstractAnemCalibr) = cal.T

"Operating calibration resistance"
resistance(cal::AbstractAnemCalibr) = cal.Rw

"Calibration fluid"
fluid(cal::AbstractAnemCalibr) = cal.fluid

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

"Tmeperature dependent resistor used for the calibration"
resistor(cal::AbstractAnemCalibr) = cal.R

struct TempCalibr{X,RT<:AbstractResistor,Fluid,Fit} <: AbstractAnemCalibr
    "Resistance element"
    R::RT
    "Calibration operating resistance"
    Rw::X
    "Calibration operating resistance temperature"
    Tw::X
    "Fluid calibration temperature"
    T::X
    "Fluid calibration pressure"
    P::X
    "Calibration fluid"
    fluid::Fluid
    "Calibration data"
    data::Matrix{X}
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
                    Ec::AbstractVector,
                    Uc::AbstractVector,
                    Tc, Pc, Rwc, makefitfun;
                    Rw=nothing, T=nothing, P=nothing,
                    fluid=AIR) where {RT<:AbstractResistor}
    
    Tm, Pm, caltab = makecaltable(Ec, Uc, Tc, Pc)

    # Reference values
    if isnothing(T)
        T = Tm
    end
    if isnothing(P)
        P = Pm
    end
    if isnothing(Rw)
        Rw = mean(Rwc)
    end

    Twc = temperature.(R, Rwc)
    
    Tw = temperature(R, Rw)

    # Correct voltages to conditions Rw, T, P
    E = Ec .* sqrt.( Rw./Rwc .* (Tw - T) ./ (Twc .- Tc)  )
    
    fit = makefitfun(E, Uc)

    
    return TempCalibr(R, Rw, Tw, T, P, fluid, caltab, fit)
    
end

function hwcorrect(cal::TempCalibr; T=caltemp(cal), P=pressure(cal),
                   fluid=fluid(cal),Rw=resistance(cal))
    
    Tw = temperature(cal.R, Rw)
    ef = sqrt( cal.Rw/Rw * (cal.Tw - cal.T) / (Tw - T) )
    return ef, one(ef)
end

"""
`velocity(cal, E; T, P, fluid, Rw)`

Calculate the velocity given anemometer output `E` and flow conditions `T`, `P`, with fluid `fluid` and operating at resistance `Rw`.

"""    
function velocity(cal::TempCalibr, E; T=caltemp(cal), P=pressure(cal),
                  fluid=fluid(cal),Rw=resistance(cal))
    fcorr,uf = hwcorrect(cal; T=T, P=P, fluid=fluid, Rw=Rw)
    return uf*cal.fit(fcorr * E)
end

"""
`velocity!(U, cal, E; T, P, fluid, Rw)`

Calculate the velocity given anemometer output `E` and flow conditions `T`, `P`, with fluid `fluid` and operating at resistance `Rw` and store it in preallocated `U`.

"""    
function velocity!(U::AbstractArray, cal::TempCalibr, E::AbstractArray;
                   T=caltemp(cal), P=pressure(cal),
                   fluid=fluid(cal),Rw=resistance(cal))
    @assert size(U) == size(E)
    fcorr,uf = hwcorrect(cal; T=T, P=P, fluid=fluid, Rw=Rw)
    for (i,e) in enumerate(E)
        U[i] = uf*cal.fit(e * fcorr)
    end
    return U
        
end

velocity(cal::TempCalibr, E::AbstractArray;
         T=temperature(cal), P=pressure(cal),
         fluid=fluid(cal),Rw=resistance(cal)) = velocity!(similar(E), cal, E; T=T,
                                                          P=P, fluid=fluid, Rw=Rw)

struct HWCalibr{X<:AbstractFloat,RT<:AbstractResistor,
                         Fluid,Fit} <: AbstractAnemCalibr
    "Resistance element"
    R::RT
    "Calibration operating resistance"
    Rw::X
    "Calibration operating resistance temperature"
    Tw::X
    "Fluid calibration temperature"
    T::X
    "Fluid calibration pressure"
    P::X
    "Calibration fluid"
    fluid::Fluid
    "Calibration data"
    data::Matrix{X}
    "Calibration curve fit"
    fit::Fit
    "Prandtl number exponent in Nusselt"
    n::X
    "Mean surface temperature factor"
    theta::X
    "Calibration kinematic viscosity"
    nu::X
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
                  Ec::AbstractVector,
                  Uc::AbstractVector,
                  Tc, Pc, Rwc, makefitfun; n=0.3, theta=0.5,
                  Rw=nothing, T=nothing, P=nothing,
                  fluid=AIR) where {RT<:AbstractResistor}
    
    Tm, Pm, caltab = makecaltable(Ec, Uc, Tc, Pc)
    
    # Reference values
    if isnothing(T)
        T = Tm
    end
    if isnothing(P)
        P = Pm
    end
    if isnothing(Rw)
        Rw = mean(Rwc)
    end
    
    Twc = temperature.(R, Rwc)
    
    Tw = temperature(R, Rw)
    ΔTc = Twc .- Tc
    ΔT  = Tw - T
    
    # Film temperature
    Tfc = Tc .+ theta .* ΔTc  # Individual calibration points
    Tf  = T   + theta  * ΔT  # Calibration reference conditions
    phic   = heatcond.(fluid, Tfc,  Pc) .* prandtl.(fluid, Tfc, Pc) .^ n

    # Reynolds number
    nuc = kinvisc.(fluid, Tfc, Pc)
    nu  = kinvisc(fluid, Tf, P)
    Rec = Uc ./ nuc

    # E²/ϕΔT
    Ex = Ec .* Ec ./ (phic .* Rw .* (Twc .- Tc))

    fit = makefitfun(Ex, Rec)
    return HWCalibr(R, Rw, Tw, T, P, fluid, caltab, fit, n, theta, nu)
    
end


function velocity(cal::HWCalibr, E; T=caltemp(cal), P=pressure(cal),
                  fluid=fluid(cal),Rw=resistance(cal))
    
    Tw = temperature(cal.R, Rw)

    ΔT = Tw - T
    Tf = T + cal.theta * ΔT
    
    k = heatcond(fluid, Tf, P)
    Pr = prandtl(fluid, Tf, P)
    phi = k*Pr^cal.n

    nu = kinvisc(fluid, Tf, P)
    
    return nu * cal.fit(E*E/(phi * Rw * ΔT))
end



function velocity!(U::AbstractArray, cal::HWCalibr, E::AbstractArray;
                   T=caltemp(cal), P=pressure(cal),
                   fluid=fluid(cal),Rw=resistance(cal))

    @assert size(E) == size(U)
    
    Tw = temperature(cal.R, Rw)
    ΔT = Tw - T
    Tf = T + cal.theta * ΔT
    k = heatcond(fluid, Tf, P)
    Pr = prandtl(fluid, Tf, P)
    den = k*Pr^cal.n * Rw * (Tw - T)
    
    nu = kinvisc(fluid, Tf, P)

    for (i,e) in enumerate(E)
        U[i] = nu * cal.fit(e*e/den)
    end
    return U
end


velocity(cal::HWCalibr, E::AbstractArray;
         T=caltemp(cal), P=pressure(cal),
         fluid=fluid(cal),Rw=resistance(cal)) = velocity!(similar(E), cal, E; T=T,
                                                             P=P, fluid=fluid, Rw=Rw)

