

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
    "Voltage offset"
    offset::U
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
offset(w::AbstractCTA) = w.offset

fluid(w::CTASensor) = fluid(w.corr)
pressure(w::CTASensor) = pressure(w.corr)
caltemp(w::CTASensor) = reftemp(w.corr)
kinvisc(w::CTASensor) = kinvisc(w.corr)
   
function correct(w::CTASensor, E;
                 T=caltemp(w), P=pressure(w),
                 fluid=fluid(w), Rw=resistance(w))
    if Rw == resistance(w)
        Tw = temperature(w)
    else
        Tw = temperature(resistor(w), Rw)
    end
    g = gain(w)
    offset = w.offset
    E1 = E/g + offset
    return correct(E1, w.corr, T, P, fluid, Rw, Tw)
    
end


#velocity(w::CTASensor, E) = w.fit(E)

function velocity(w::CTASensor, E::Real;
                 T=caltemp(w), P=pressure(w),
                  fluid=fluid(w), Rw=resistance(w))
    fc = correct(w, E; T=T, P=P, fluid=fluid, Rw=Rw)
    return velocity(w, E, fc)
end

function velocity(w::CTASensor, E::Real, fc::CorrFactor)
    Ec = fc(E, gain(w), offset(w))
    Uc = w.fit(Ec)
    return (fc.nu / kinvisc(w.corr)) * Uc
end

function velocity!(U::AbstractVector, w::CTASensor, E::AbstractVector, fc::CorrFactor)
    @assert length(U) == length(E)
    
    g = gain(w)
    o = offset(w)

    rν = fc.nu/kinvisc(w.corr)
    
    for (i,e) in enumerate(E)
        ec = fc(e, g, o)
        uc = w.fit(ec)
        U[i] = rν * uc
    end
    return U
end
velocity(w::CTASensor, E::AbstractVector, fc::CorrFactor) =
    velocity(zeros(length(E)), w, W, fc)

function velocity!(U::AbstractVector, w::CTASensor, E::AbstractVector;
                   T=caltemp(w), P=pressure(w),
                   fluid=fluid(w), Rw=resistance(w))
    fc = correct(w, E; T=T, P=P, fluid=fluid, Rw=Rw)
    return velocity!(U, w, E, fc)
end

function velocity!(w::CTASensor, E::AbstractVector;
                   T=caltemp(w), P=pressure(w),
                   fluid=fluid(w), Rw=resistance(w))
    fc = correct(w, E; T=T, P=P, fluid=fluid, Rw=Rw)
    return velocity!(zeros(length(E)), w, E, fc)
end


(w::CTASensor)(E; kw...) = velocity(w, E; kw...)
(w::CTASensor)(E, fc::CorrFactor) = velocity(w, E, fc)


