

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
struct CCASensor{Correct,U,RT,Fit} <: AbstractCCA
    "Temperature dependent resistor"
    R::RT
    "Operating electrical current"
    I::U
    "Voltage output gain"
    gain::U
    "Correction model"
    corr::Correct
    "Calibration curve"
    fit::Fit
end

Base.broadcastable(sensor::CCASensor) = Ref(sensor)

resistor(w::AbstractCCA) = w.R

"Operating current of the CTA"
current(w::AbstractCCA) = w.Rw

reftemp(w::CCASensor) = reftemp(w.R)
refresist(w::CCASensor) = refresist(w.R)

gain(w::AbstractCCA) = w.gain

fluid(w::CCASensor) = fluid(w.corr)
pressure(w::CCASensor) = pressure(w.corr)
caltemp(w::CCASensor) = reftemp(w.corr)


function correct(w::CCASensor, E;
                 T=caltemp(w), P=pressure(w),
                 fluid=fluid(w), I=current(w))
    g = gain(w)
    Ew = E / g # Voltage passing through the resistor
    Rw = Ew / I # Resistance of the sensor
    Tw = temperature(resistor(w), Rw)
    
    return correct(Ew, w.corr, T, P, fluid, Rw, Tw)*g
        
end

function velocity(w::CCASensor, E;
                  T=caltemp(w), P=pressure(w),
                  fluid=fluid(w), I=current(w))
    g = gain(w)
    Ew = E / g # Voltage passing through the resistor
    Rw = Ew / I # Resistance of the sensor
    Tw = temperature(resistor(w), Rw)
    
    Ec = correct(Ew, w.corr, T, P, fluid, Rw, Tw)*g
    Uc = w.fit(Ec)
    ν_cal = kinvisc(w.corr)
    ν     = kinvisc(fluid, (T+Tw)/2, P)
    return ν / ν_cal * Uc
    
end
    
