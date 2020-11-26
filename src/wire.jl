
abstract type AbstractThermalAnemometer end

"""
    CTASensor(R, a)

A structure to manage constant temperature anemometer sensors (CTA)

 * `R` Resistance element
 * `a` overheat ratio

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
    Tw::Float64
    "Reference temperature"
end

"Resistor element of the CTA sensor"
resistor(w::CTASensor) = w.R

"Operating temperature of the CTA"
optemperature(w::CTASensor) = w.Tw

"Overheat ratio of the CTA sensor"
overheat_ratio(w::CTASensor) = resistance(w.R, w.Tw)/resistance(w.R) - 1.0


"Over temperature"
overtemp(w::CTASensor) = optemperature(w) - temperature(w.R)


                                          
                                          
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
end

resistor(w::CCASensor) = w.R
"Current in A flowing through the CCA"
current(w::CCASensor) = w.I
