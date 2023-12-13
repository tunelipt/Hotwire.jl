

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
struct CTASensor{Calibr,U,RT} <: AbstractCTA
    "Temperature dependent resistor"
    R::RT
    "Operating resistance of the sensor"
    Rw::U
    "Operating temperature of the sensor"
    Tw::U
    "Voltage output gain"
    gain::U
    "Calibration curve U = cal(E)"
    cal::Calibr
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)

#CTASensor(R::RT, Rw, gain, cal, corr,fluid=AIR) where {RT<:AbstractResistor} =
#
CTASensor(R, Rw, temperature(R,Rw), gain, cal, corr)
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

calibr(w::AbstsractThermalAnemometer) = w.cal
fluid(c::AbstractThermalAnemometer) = c.fluid
pressure(w::AbstractThermalAnemometer) = pressure(calibr(w))

function correction(w::CTASensor{AC}, E;
                    T=reftemp(w.cal), Rw=resistance(w),
                    fluid=fluid(w.cal), P=pressure(w.cal)) where {AC}
    g = gain(w)
    # Voltage accross resistive element
    E1 = E / g
    
    R = resistor(w)
    Tw = temperature(R, Rw)

    # Operating conditions
    tc = TempCorrect(T, Rw, Tw)
    # Model at operating conditions
    correction(calibt(w), E, tc, fluid, P) * g
end

function velocity(w::CTASensor, E;
                  T=reftemp(w.cal), Rw=resistance(w),
                  fluid=fluid(w.cal), P=pressure(w.cal))
    Ec = correction(w, E; T=T, Rw=Rw, fluid=fluid, P=P)
end



velocity(w::CTASensor, E::Real) = calibr(w).fit(E)

(w::CTASensor)(E) = velocity(w, E)
(w::CTASensor)(E; T=reftemp(w.cal), Rw=resistance(w),
               fluid=fluid(w.cal), P=pressure(w.cal)) = velocity(w, E,
                                                                 T=T,Rw=Rw,
                                                                 fluid=fluid, P=P)


        
                  
                  
   
