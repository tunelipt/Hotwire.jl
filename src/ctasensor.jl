

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
struct CTASensor{Correct,Calibr,T,RT} <: AbstractCTA
    "Temperature dependent resistor"
    R::RT
    "Operating resistance of the sensor"
    Rw::T
    "Operating temperature of the sensor"
    Tw::T
    "Reference temperature of the sensor"
    Tr::T
    "Voltage output gain"
    gain::T
    "Calibration curve"
    cal::Calibr
    "Thermal model of the sensor"
    corr::Correct
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

resistor(w::AbstractCTA) = w.RT

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

gain(w::AbstractACTA) = w.gain

calibr(w::CTASensor) = w.cal

function velocity(w::CTASensor{AC}, E, meas::AC) where {AC<:AbstractAnemCorrect}
    cal = calibr(w)
    g = gain(w)
    Ec = anemcorrect(E/g, w.corr, meas)
    return  g*cal(Ec) * (kinvisc(meas) /  kinvisc(w.corr))
end

(w::CTASensor{AC})(E, meas::AC) where {AC<:AbstractAnemCorrect} = velocity(w,E,meas)


function velocity(w::CTASensor, E;
                  T=reftemp(w.corr), P=101325.0,
                  x=fluid(w.corr),R=resistance(w.corr))
    
    g = gain(w)
    if R != resistance(w.corr)
        Tw = temperature(R)
    else
        Tw = temperature(w.corr)
    end
    Ec = anemcorrect(E/g, w.corr, T, P, x, R, Tw)
    return  g*cal(Ec) * (kinvisc(meas) /  kinvisc(w.corr))
end

    
        
                  
                  
   
