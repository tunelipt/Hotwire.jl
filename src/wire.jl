
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
    "Output gain"
    gain::Float64
end

CTASensor(R::ResType, Rw, gain=1.0) where {ResType <: AbstractResistor} = CTASensor(R,Rw,gain)
"Resistor element of the CTA sensor"
resistor(w::CTASensor) = w.R

"Operating temperature of the CTA"
optemperature(w::CTASensor) = temperature(w.R, w.Rw)

"Overheat ratio of the CTA sensor"
overheat_ratio(w::CTASensor) = w.Rw/resistance(w.R) - 1.0



"Over temperature"
overtemp(w::CTASensor) = optemperature(w) - temperature(w.R)


"Gain of anemometer output"
gain(w::AbstractThermalAnemometer) = w.gain

                                          
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
end
CCASensor(R::ResType, I, gain=1.0) where {ResType <: AbstractResistor} = CCASensor(R,I,gain)

resistor(w::CCASensor) = w.R
"Current in A flowing through the CCA"
current(w::CCASensor) = w.I
