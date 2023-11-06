```@meta
CurrentModule = Hotwire
```

# Hotwire: a [Julia](https://julialang.org) package for thermal anemometry

Thermal anemometry is a class of instruments used to measure fluid velocities based on convective heat transfer from a heated source to a fluid. It isa very flexible methodology to measure fluid velocities but it does have its quirks.

The [Hotwire.jl](https://github.com/tunelipt/Hotwire.jl) package deals sensors whose electrical resistance varies with resistor temperature due to self heating.

The two most common modes of operation are:

 * Constant temperature anemometer (CTA)
 * Constant current anemometer (CCA)

This package deals with both of these types of sensors even though emphasis is on  constant temperature sensors because of their wider use.

The package implements two types of resistors

 * Linear resistors
 * NTC Thermistors

But the interface tries to be a generic as possible such that other types of resistors can be used if the correct interface is used.

The name Hotwire comes from one of the most important and earliest applications of thermal anemometry: the [hotwire anemometer](https://en.wikipedia.org/wiki/Anemometer#Hot-wire_anemometers). It was the first sensor to be capable of measuring turbulent velocity fluctuations and is still today the most used type of sensor when measuring turbulent flows.

# Basic principles

The thermal anemometer uses the convective heat transfer from a heated resistor to estimate the velocity. The resistor is heated due to an electric current passing through the resistor. Given a heated body, the heat loss to a fluid is modeled by [Newton's law of cooling](https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling):

``\dot{Q} = h \cdot A \cdot \left( T_w - T_a \right)``





Documentation for [Hotwire](https://github.com/tunelipt/Hotwire.jl).

```@index
```

```@autodocs
Modules = [Hotwire]
```
