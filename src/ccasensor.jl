

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
struct CTASensor{Correct,Calibr,U,RT} <: AbstractCCA
    "Temperature dependent resistor"
    R::RT
    "Operating electrical current"
    I::U
    "Voltage output gain"
    gain::U
    "Calibration curve"
    cal::Calibr
    "Thermal model of the sensor"
    corr::Correct
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

CTASensor(R, Rw, gain, cal, corr) = CTASensor
resistor(w::AbstractCCA) = w.R

"Operating current of the CTA"
current(w::AbstractCCA) = w.Rw

reftemp(w::CCASensor) = reftemp(w.R)
refresist(w::CCASensor) = refresist(w.R)


gain(w::AbstractCCA) = w.gain

calibr(w::CCASensor) = w.cal

function velocity(w::CCASensor{AC}, E,
                  meas::AC, I=current(w)) where {AC<:AbstractAnemCorrect}
    cal = calibr(w)
    g = gain(w)
    E1 = E / g
    Rw = E1 / w.I
    Tw = temperature(w.R, Rw)
    # Something funny. I still need to implement
    Ec = anemcorrect(E/g, w.corr, meas)
    return  g*cal(Ec) * (kinvisc(meas) /  kinvisc(w.corr))
end

(w::CCASensor{AC})(E, meas::AC) where {AC<:AbstractAnemCorrect} = velocity(w,E,meas)


function velocity(w::CCASensor{Correct}, E;
                  T=reftemp(w.corr), P=101325.0,
                  x=fluid(w.corr),R=resistance(w.corr)) where {Correct}

    if R != resistance(w.corr)
        Tw = temperature(R)
    else
        Tw = temperature(w.corr)
    end
    meas = Correct(w.corr, T, Rw, Tw)
    return velocity(w, E, meas, current(w))
end

    
(w::CCASensor{Correct})(E;
                  T=reftemp(w.corr), P=101325.0,
                  x=fluid(w.corr),R=resistance(w.corr)) where {Correct}        
                  
                  
   
