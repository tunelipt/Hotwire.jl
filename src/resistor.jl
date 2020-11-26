

abstract type AbstractResistor end

"""
    r = Resistor(R0=1e3, α=0.0, T₀=20.0)

A resistor with resistance varying linearly with temperature:

```math
R = R₀ (1 + α (T - T₀))
```

**IMPORTANT** In the constructor, coefficient `α` is a *percentage*
 linear. 

# Exemples

```jldoctest
julia> R = Resistor(1e3, 0.01)
THW.Resistor(1000.0, 0.0001, 20.0)

julia> resistance(R)
1000.0

julia> resistance(R, 120)
1010.0

julia> R()
1000.0

julia> R(120)
1010.0

julia> temperature(R)
20.0

julia> temperature(R, 1010)
120.00000000000009

```
"""
struct Resistor <: AbstractResistor
    "Resistance at reference temperature"
    R0::Float64
    "Linear resistance coefficient"
    α::Float64
    "Reference temperature"
    T0::Float64
    Resistor(R0=1e3, α=0.0, T₀=20.0) = new(R0, α/100, T₀)
end




"""
    R = Thermistor(R₀, B, T₀)

Models a NTC thermistor. In the constructor, `R0` is the resistance of the thermistor at the reference temperature  `T0`. For convenience, the temperature should be expressed in °C, but it is stored in K. `B` is a coefficient that characterizes the temperature dependence of the thermistor. In the model used her, the resistance varies with temperature according to the following equation:

```math
R = R₀ \\exp B \\left(\\frac{1}{T₀} - \\frac{1}{T}\\right)
```

where
 * R₀ is the reference resistance in Ω
 * T₀ is the reference temperature in K, usually the adopted value is   293.15 K
 * T is the temperature in K, where the value of resistance should be calculated
 * B is an empirical coefficient in K that characterizes the temperature dependency of the resistance

    When using the struct `Thermistor`, beware: in the constructor, the unit of T0 is °C NOT K.

# Examples
```jldoctest
julia> R = Thermistor(5e3, 3200, 20) # Create an object `Thermistor`
Thermistor(5000.0, 3200.0, 293.15)

julia> resistance(R) # Resistance at the reference temperature
5000.0

julia> resistance(R, 25) # Resistance at 25°C
4163.587774917559

julia> R() # The same as resistance(R)
5000.0

julia> R(25) # The same as resistance(R, 25)
4163.587774917559

julia> temperature(R)
20.0

julia> temperature(R, 4163.588)
24.999998498264233

```

"""
struct Thermistor <: AbstractResistor
    "Reference resistance in Ω at temperature `T0`"
    R0::Float64
    "Thermistor's B coefficient in K "
    B::Float64
    "Reference temperature in K"
    T0::Float64
    Thermistor(R0=5e3, B=0.0, T0=25.0) = new(R0, B, T0+273.15)
end



"""
    resistance(th::AbstractResistor)
    resistance(th::AbstractResistor, Tc)

Calculates the resistance for `AbstractResistor`objects` at temperature `Tc`.
If argumento `Tc` is not provided, the function returns the reference resistance `R0`.

See [`Thermistor`](@ref) to see more details of the  type  `Thermistor`. For more details on 
type `Resistor`, see [`Thermistor`](@ref).

"""
resistance(r::AbstractResistor) = r.R0
resistance(r::Resistor, T) = r.R0 * (1.0 + r.α * (T - r.T0))
resistance(th::Thermistor, T) = th.R0 * exp( -th.B * (1/th.T0 - 1/(T+273.15) ) )

temperature(r::AbstractResistor) = r.T0
temperature(r::Resistor, R) = 1/r.α * (R/r.R0 - 1.0) + r.T0
temperature(th::Thermistor, R) = 1/( 1/th.T0 + 1/th.B * log(R/th.R0) ) - 273.15

(th::Thermistor)(T) = resistance(th, T)
(th::Thermistor)() = th.R0
(r::Resistor)(T) = resistance(r, T)
(r::Resistor)() = resistance(r)

    
