

"""
    CCASensor(R, I)

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
struct CCASensor{U,RT} <: AbstractCCA
    "Temperature dependent resistor"
    R::RT
    "Operating electrical current"
    I::U
    "Voltage output gain"
    gain::U
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

CCASensor(R, Rw, gain, cal, corr) = CCASensor
resistor(w::AbstractCCA) = w.R

"Operating current of the CTA"
current(w::AbstractCCA) = w.Rw

reftemp(w::CCASensor) = reftemp(w.R)
refresist(w::CCASensor) = refresist(w.R)


gain(w::AbstractCCA) = w.gain

