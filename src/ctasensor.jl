



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
    "Operating temperature of the sensor"
    Tw::U
    "Calibration and correction model"
    calibr::Calibr
end
Base.broadcastable(sensor::CTASensor) = Ref(sensor)

CTASensor(calibr) = CTASensor(calibr.R, calibr.Rw, calibr.Tw, calibr)

function CTASensor(calibrfun, R::RT, Rw, Ec, Uc, Tc, Pc, makefitfun;
                   T=nothing, P=nothing, fluid=AIR, kw...) where {RT<:AbstractResistor}
    calibr = calibrfun(R, Ec, Uc, Tc, Pc, Rw, makefitfun; kw...)
    return CTASensor(calibr)
end

resistor(w::AbstractCTA) = w.R

"Operating resistance of the CTA"
resistance(w::AbstractCTA) = w.Rw

"Operating temperature of the CTA"
temperature(w::CTASensor) = w.Tw

reftemp(w::CTASensor) = reftemp(w.R)
refresist(w::CTASensor) = refresist(w.R)

"Overheat ratio of the CTA"
overheatratio(w::AbstractCTA, T) =
    resistance(w) / resistance(resistor(w),T) - 1

overheatratio(w::AbstractCTA) = resistance(w) / refresist(w) - 1

"Temperature above reference temperature that the CTA sensor operates"
overtemp(w::AbstractCTA,T) = temperature(w) - T
overtemp(w::AbstractCTA) = temperature(w) - reftemp(w) 


fluid(w::CTASensor) = fluid(w.calibr)
pressure(w::CTASensor) = pressure(w.calibr)
caltemp(w::CTASensor) = reftemp(w.calibr)
kinvisc(w::CTASensor) = kinvisc(w.calibr)




function velocity(w::CTASensor, E::Real;
                 T=caltemp(w), P=pressure(w),
                  fluid=fluid(w), Rw=resistance(w))
    
    return velocity(w.calibr, E; T=T, P=P, fluid=fluid, Rw=Rw)
end


function velocity!(U::AbstractArray, w::CTASensor, E::AbstractArray;
                   T=caltemp(w), P=pressure(w),
                   fluid=fluid(w), Rw=resistance(w))
    
    return velocity!(U, w.calibr, E; T=T, P=P, fluid=fluid, Rw=Rw)
end

velocity(w::CTASensor, E::AbstractArray;
         T=caltemp(w), P=pressure(w),
         fluid=fluid(w), Rw=resistance(w)) = velocity!(similar(E), w.calibr, E,
                                                       T=T,P=P,fluid=fluid,Rw=Rw)


(w::CTASensor)(E; kw...) = velocity(w, E; kw...)


