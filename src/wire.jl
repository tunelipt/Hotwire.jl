
abstract type AbstractThermalAnemometer end

"""
    CTASensor(R, Rw)

A structure to manage constant temperature anemometer sensors (CTA)

 * `R` Resistance element
 * `Rw` Operating resistance
 * `gain` Output gain, ``Eo = g⋅E``

The overheat ratio is defined by the ratio between over-resistance and reference 
resistance:

    `a = (Rw - Ro) / Ro`
where `a` is the overheat ratio, `Rw` is the operating resistance of the element sensor, 
`Ro` is the reference resistance (resistance at reference temperature).

"""
struct CTASensor{ResType <: AbstractResistor} <: AbstractThermalAnemometer
    "Element whose resistance changes with temperature"
    R::ResType
    "Operating temperature of the sensor"
    Rw::Float64
    "Reference temperature"
    T₀::Float64
    "Output gain"
    gain::Float64
    
end
Base.broadcastable(sensor::CTASensor) = Ref(sensor)

CTASensor(R::ResType, Rw; gain=1.0) where {ResType <: AbstractResistor} = CTASensor(R, Rw, temperature(R), gain)

CTASensor(R::ResType, Rw, T₀; gain=1.0) where {ResType <: AbstractResistor} = CTASensor(R, Rw, T₀, gain)


Wire(R₀, Rw, T₀=293.15, α=0.36e-3, gain=1.0) = CTASensor(Resistor(R₀, α, T₀), Rw, T₀; gain=gain)

"Operating temperature of the CTA"
optemperature(w::CTASensor) = temperature(w.R, w.Rw)

"Overheat ratio of the CTA sensor"
overheat_ratio(w::CTASensor) = w.Rw/resistance(w.R) - 1.0



"Over temperature"
overtemp(w::CTASensor) = optemperature(w) - temperature(w.R)


"Resistor element of a thermal anemometer sensor"
resistor(w::AbstractThermalAnemometer) = w.R

"Gain of anemometer output"
gain(w::AbstractThermalAnemometer) = w.gain

"Reference temperature of sensor"
reftemp(w::CTASensor) = w.T₀

"""
    CCASensor(R, a)

A structure to manage constant current anemometer sensors (CCA)

 * `R` Resistance element
 * `I` current in Ampere


"""
struct CCASensor{ResType <: AbstractResistor} <: AbstractThermalAnemometer
    "Element whose resistance changes with temperature"
    R::ResType
    "Operating current in Ampere"
    I::Float64
    "Output gain"
    gain::Float64
    "Reference temperature"
    T₀::Float64
end
CCASensor(R::ResType, I, T0=reftemp(R); gain=1.0) where {ResType <: AbstractResistor} = CCASensor(R,I,gain,T0)
Base.broadcastable(sensor::CCASensor) = Ref(sensor)

"Current in A flowing through the CCA"
current(w::CCASensor) = w.I

reftemp(w::CCASensor) = w.T₀
