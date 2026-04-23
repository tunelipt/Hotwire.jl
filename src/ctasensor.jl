



"""
    CTASensor(R, Rw)

A structure to manage constant temperature anemometer sensors (CTA)

 * `R` Resistance element
 * `Rw` Operating resistance
 * `gain` Output gain, ``Eo = g⋅E``

The overheat ratio is defined by the ratio between over-resistance and reference 
resistance:

    ``a = (R_w - R₀) / R₀`
where `a` is the overheat ratio, `Rw` is the operating resistance of the element sensor, 
`Ro` is the reference resistance (resistance at reference temperature).

"""
struct CTASensor{U,RT,Calibr} <: AbstractCTA
    "Temperature dependent resistor"
    R::RT
    "Operating resistance of the sensor"
    Rw::U
    "Calibration and correction model"
    calibr::Calibr
end
Base.broadcastable(sensor::CTASensor) = Ref(sensor)




resistor(w::AbstractCTA) = w.R

"Operating resistance of the CTA"
resistance(w::AbstractCTA) = w.Rw

"Operating temperature of the CTA"
temperature(w::CTASensor) = temperature(w.R, w.Rw)

reftemp(w::CTASensor) = reftemp(w.R)
refresist(w::CTASensor) = refresist(w.R)

"Overheat ratio of the CTA"
overheatratio(w::AbstractCTA, T) =
    resistance(w) / resistance(resistor(w),T) - 1

overheatratio(w::AbstractCTA) = resistance(w) / refresist(w) - 1

"Temperature above reference temperature that the CTA sensor operates"
overtemp(w::AbstractCTA,T) = temperature(w) - T
overtemp(w::AbstractCTA) = temperature(w) - reftemp(w) 


calibration(w::AbstractCTA) = w.calibr

caltemp(w) = caltemp(calibration(w))
calpress(w) = calpress(calibration(w))
calfluid(w) = calfluid(calibration(w))


function velocity(w::CTASensor, E; T=caltemp(w), P=calpress(w),
                  fluid=calfluid(w), R=resistance(Rw))
    velocity(calibration(w), resistor(w), E, R, T, P, fluid)
end

function velocity!(U::AbstractArray, w::CTASensor, E::AbstractArray;
                   T=caltemp(w), P=calpress(w),
                   fluid=calfluid(w), R=resistance(Rw))
    @assert size(U) == size(E)
    param, nu = hwcorrect(calibration(w), resistor(w), R, T, P, fluid)
    map!(e->nu * calibration(w).fit(e*e/param), U, E)
end

function velocity(w::CTASensor, E::AbstractArray; T=caltemp(w), P=calpress(w),
                  fluid=calfluid(w), R=resistance(Rw))
    velocity!(similar(W), w, E; T=T, P=P, fluid=fluid, R=R)
end


(w::CTASensor)(E; kw...) = velocity(w, E; kw...)


