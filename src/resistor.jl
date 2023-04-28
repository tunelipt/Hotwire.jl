
"""
`AbstractResistor`

Base abstract class for defining temperature dependent resistors.

An `AbstractResistor` should implement the following methods:

 * [`resistance`](@ref) calculates the resistance at a specific temperature.
 * [`temperature`](@ref) calculates the temperature where the resistor attains the specified resistance.


"""
abstract type AbstractResistor{T} end

"""
    r = Resistor(R0=1e3, α=0.0, T₀=293.15)

A resistor with resistance varying linearly with temperature:

```math
R = R₀ (1 + α (T - T₀))
```

**IMPORTANT** In the constructor, coefficient `α` is a *fraction*, not percentage.
 linear. 

# Examples

```jldoctest
julia> R = Resistor(1e3, 0.01)
THW.Resistor(1000.0, 0.0001, 293.15)

julia> resistance(R)
1000.0

julia> resistance(R, 120+273.15)
1010.0

julia> R()
1000.0

julia> R(120+273.15)
1010.0

julia> temperature(R)
293.15

julia> temperature(R, 1010)
120.00000000000009

```
"""
struct Resistor{T<:Number} <: AbstractResistor{T}
    "Resistance at reference temperature (Ω)"
    R₀::T
    "Linear resistance coefficient"
    α::T
    "Reference temperature (K)"
    T₀::T
end
Base.broadcastable(R::Resistor) = Ref(R)

Resistor(R₀::T, α::U, T₀::V) where {T,U,V} =
    Resistor(promote(R₀,α,T₀)...)

Resistor(;R0,a=0.4e-2,T0=20.0) = Resistor(R0,a,T0)


"""
    R = Thermistor(R₀, B, T₀)

Models a NTC thermistor. In the constructor, `R0` is the resistance of the thermistor at the reference temperature  `T₀`. For convenience, the temperature should be expressed in K. `B` is a coefficient that characterizes the temperature dependence of the thermistor. In the model used her, the resistance varies with temperature according to the following equation:

```math
R = R₀ \\exp B \\left(\\frac{1}{T₀} - \\frac{1}{T}\\right)
```

where
 * R₀ is the reference resistance in Ω
 * T₀ is the reference temperature in K, usually the adopted value is   293.15 K
 * T is the temperature in K, where the value of resistance should be calculated
 * B is an empirical coefficient in K that characterizes the temperature dependency of the resistance

    When using the struct `Thermistor`, beware: in the constructor, the unit of T₀ is °C NOT K.

# Examples
```jldoctest
julia> R = Thermistor(5e3, 3200, 293.15) # Create an object `Thermistor`
Thermistor(5000.0, 3200.0, 293.15)

julia> resistance(R) # Resistance at the reference temperature
5000.0

julia> resistance(R, 298.15) # Resistance at 25°C
4163.587774917559

julia> R() # The same as resistance(R)
5000.0

julia> R(298.15) # The same as resistance(R, 298.15)
4163.587774917559

julia> temperature(R)
293.15

julia> temperature(R, 4163.588)
298.15

```

"""
struct Thermistor{T<:Number} <: AbstractResistor{T}
    "Reference resistance in Ω at temperature `T₀`"
    R₀::T
    "Thermistor's B coefficient in K "
    B::T
    "Reference temperature in K"
    T₀::T
end
Base.broadcastable(R::Thermistor) = Ref(R)
Thermistor(R₀::T, B::U, T₀::V) where {T,U,V} =
    Thermistor(promote(R₀, B, T₀)...)

Thermistor(;R0=5e3,B=3950.0,T0=298.15) = Thermistor(R0,B,T0)



"""
    `resistance(th::AbstractResistor)`
    `resistance(th::AbstractResistor, T)`

Calculates the resistance for `AbstractResistor`objects` at temperature `T`.
If argument `T` is not provided, the function returns the reference resistance `R0`.

See [`Thermistor`](@ref) to see more details of the  type  `Thermistor`.
For more details on type `Resistor`, see [`Thermistor`](@ref).

See also [`temperature`](@ref), [`reftemp`](@ref), [`refresist`](@ref)

"""
resistance(r::AbstractResistor) = r.R₀
resistance(r::Resistor, T) = r.R₀ * (1 + r.α * (T - r.T₀))
resistance(th::Thermistor, T) = th.R₀ * exp( -th.B * (1/th.T₀ - 1/T ) )

"""
    `reftemp(R)`

Returns the reference temperature of an [`AbstractResistor`](@ref) object.
See [`refresist`](@ref) for returning the corresponding reference resistance.

Method [`temperature`](@ref) returns the temperature for a given resistance value
and methods [`resistance`](@ref) returns the resistance at a given temperature.
"""
reftemp(r::AbstractResistor) = r.T₀

"""
    `refresist(R)`

Returns the reference resistance of an [`AbstractResistor`](@ref) object.
See [`reftemp`](@ref) for returning corresponding the reference temperature.
"""
refresist(r::AbstractResistor) = r.R₀

"""
    `temperature(R)`
    `temperature(R, Rx)`

Returns the temperature at which an [`AbstractResistor`](@ref) reaches a resistance
`Rx`. If `Rx` is not provided, the method returns the reference temperature.

See [`Thermistor`](@ref) to see more details of the  type  `Thermistor`.
For more details on type `Resistor`, see [`Thermistor`](@ref).

See also [`resistance`](@ref), [`reftemp`](@ref), [`refresist`](@ref)
"""
temperature(r::AbstractResistor) = r.T₀
temperature(r::Resistor, R) = 1/r.α * (R/r.R₀ - 1) + r.T₀
temperature(th::Thermistor, R) = 1/( 1/th.T₀ + 1/th.B * log(R/th.R₀) ) 
temperature(th::Thermistor) = th.T₀ 

(th::Thermistor)(T) = resistance(th, T)
(th::Thermistor)() = th.R₀
(r::Resistor)(T) = resistance(r, T)
(r::Resistor)() = resistance(r)

    
