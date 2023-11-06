

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
struct CTASensor{T,U,RT,Calibr,Correct} <: AbstractCTA
    R::RT
    Rw::T
    gain::U
    cal::Calibr
    corr::Correct
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

resistor(w::AbstractCTA) = w.R
resistance(w::AbstractCTA) = w.Rw
temperature(w::AbstractCTA) = temperature(sensor(w), resistance(w))
overheatratio(w::AbstractCTA, T) = resistance(w) / resistance(sensor(w),T) - 1
overheatratio(w::AbstractCTA) = resistance(w) / resistance(sensor(w)) - 1
overtemp(w::AbstractCTA,T) = temperature(sensor(w),resistance(w)) - T
overtemp(w::AbstractCTA) = temperature(sensor(w),resistance(w)) - temperature(sensor(w))



