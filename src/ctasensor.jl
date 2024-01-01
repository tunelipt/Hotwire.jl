

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
struct CTASensor{Correct,U,RT,Fit} <: AbstractCTA
    "Temperature dependent resistor"
    R::RT
    "Operating resistance of the sensor"
    Rw::U
    "Operating temperature of the sensor"
    Tw::U
    "Voltage output gain"
    gain::U
    "Correction model"
    corr::Correct
    "Calibration curve"
    fit::Fit
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

resistor(w::AbstractCTA) = w.R

"Operating resistance of the CTA"
resistance(w::AbstractCTA) = w.Rw

"Operating temperature of the CTA"
temperature(w::CTASensor) = w.Tw

reftemp(w::CTASensor) = reftemp(w.R)
refresist(w::CTASensor) = refresist(w.R)

"Overheat ratio of the CTA"
overheatratio(w::AbstractCTA, T) = resistance(w) / resistance(resistor(w),T) - 1
overheatratio(w::AbstractCTA) = resistance(w) / refresist(w) - 1

"Temperature above reference temperature that the CTA sensor operates"
overtemp(w::AbstractCTA,T) = temperature(w) - T
overtemp(w::AbstractCTA) = temperature(w) - reftemp(w) 

gain(w::AbstractCTA) = w.gain
fluid(w::CTASensor) = fluid(w.corr)
pressure(w::CTASensor) = pressure(w.corr)
caltemp(w::CTASensor) = reftemp(w.corr)

   
function correct(w::CTASensor, E;
                 T=caltemp(w), P=pressure(w),
                 fluid=fluid(w), Rw=resistance(w))
    if Rw == resistance(w)
        Tw = temperature(w)
    else
        Tw = temperature(resistor(w), Rw)
    end
    g = gain(w)
    return correct(E/g, w.corr, T, P, fluid, Rw, Tw)
    
end



function velocity(w::CTASensor, E;
                 T=caltemp(w), P=pressure(w),
                 fluid=fluid(w), Rw=resistance(w))
    if Rw == resistance(w)
        Tw = temperature(w)
    else
        Tw = temperature(resistor(w))
    end
    g = gain(w)
    (fc,ν) = correct(E/g, w.corr, T, P, fluid, Rw, Tw)
    
    ν_cal = kinvisc(w.corr)
    Uc = w.fit(E*fc)
    return ν / ν_cal * Uc
end

function velocity(w::CTASensor, E, (fc,ν))
    ν_cal = kinvisc(w.corr)
    Uc = w.fit(E*fc)
    return ν/ν_cal * Uc
end

(w::CTASensor)(args...; kw...) = velocity(w, args...; kw...)

#T=caltemp(w), P=pressure(w),
#                  fluid=fluid(w), Rw=resistance(w)) =
#                      velocity(w, E; T=T, P=P, fluid=fluid,
#                               Rw=Rw)

