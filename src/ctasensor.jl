

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
struct CTASensor{Correct,Calibr,T,U,RT} <: AbstractCTA
    "Temperature dependent resistor"
    R::RT
    "Operating resistance of the sensor"
    Rw::T
    "Voltage output gain"
    gain::U
    "Calibration curve"
    cal::Calibr
    "Fluid correctio routine"
    corr::Correct
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

"Operating resistance of the CTA"
resistance(w::AbstractCTA) = w.Rw

"Operating temperature of the CTA"
temperature(w::AbstractCTA) = temperature(sensor(w), resistance(w))

"Overheat ratio of the CTA"
overheatratio(w::AbstractCTA, T) = resistance(w) / resistance(sensor(w),T) - 1
overheatratio(w::AbstractCTA) = resistance(w) / resistance(sensor(w)) - 1

"Temperature above reference temperature that the CTA sensor operates"
overtemp(w::AbstractCTA,T) = temperature(sensor(w),resistance(w)) - T
overtemp(w::AbstractCTA) = temperature(sensor(w),resistance(w)) - temperature(sensor(w))

correction(w::CTASensor) = w.corr



   
