

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
struct CTASensor{T,RT,Calibr} <: AbstractCTA
    "Temperature dependent resistor"
    R::RT
    "Operating resistance of the sensor"
    Rw::T
    "Operating temperature of the sensor"
    Tw::T
    "Reference temperature of the sensor"
    Tr::T
    "Reference resistance of the sensor"
    Rr::T
    "Voltage output gain"
    gain::T
    "Calibration curve"
    cal::Calibr
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

resistor(w::AbstractCTA) = w.RT

"Operating resistance of the CTA"
resistance(w::AbstractCTA) = w.Rw

"Operating temperature of the CTA"
temperature(w::CTASensor) = w.Tw

reftemp(w::CTASensor) = w.Tr
refresist(w::CTASensor) = w.Rr

"Overheat ratio of the CTA"
overheatratio(w::AbstractCTA, T) = resistance(w) / resistance(resistor(w),T) - 1
overheatratio(w::AbstractCTA) = resistance(w) / refresist(w) - 1

"Temperature above reference temperature that the CTA sensor operates"
overtemp(w::AbstractCTA,T) = temperature(w) - T
overtemp(w::AbstractCTA) = temperature(w) - reftemp(w) 

gain(w::AbstractACTA) = w.gain

calibr(w::CTASensor) = w.cal



   
