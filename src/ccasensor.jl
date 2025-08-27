

"""
    CCASensor(R, I)

A structure to manage constant temperature anemometer sensors (CTA)

 * `R` Resistance element
 * `I` Operating constant current

The overheat ratio is defined by the ratio between over-resistance and reference 
resistance:

    `a = (Rw - Ro) / Ro`
where `a` is the overheat ratio, `Rw` is the operating resistance of the element sensor, 
`Ro` is the reference resistance (resistance at reference temperature).

"""
struct CCASensor{Calibr,U,RT,Fit} <: AbstractCCA
    "Temperature dependent resistor"
    R::RT
    "Operating electrical current"
    I::U
    "Voltage output gain"
    gain::U
    "Calibration and correction model"
    calibr::Calibr
end
Base.broadcastable(sensor::CCASensor) = Ref(sensor)

function CCASensor(model::Calibr, R::RT, I, gain, E::AbstractArray, U::AbstractArray,
                   T, P, makefitfun;
                   fluid=AIR, params...) where {RT<:AbstractResistor,
                                                Calibr<:AbstractAnemCalibr}

    Ew = E ./ gain  # Voltage on the sensor element itself
    Rwc = Ew ./  I  # Sensor resistance
    calibr = model(R, Ew, U, T, P, Rwc, makefitfun; fluid=fluid, params...)
    # Now create the actual object
    return CCASensor(R, I, gain, calibr)

end
resistor(w::AbstractCCA) = w.R

"Operating current of the CTA"
current(w::AbstractCCA) = w.Rw

reftemp(w::CCASensor) = reftemp(w.R)
refresist(w::CCASensor) = refresist(w.R)

gain(w::AbstractCCA) = w.gain

fluid(w::CCASensor) = fluid(w.calibr)
pressure(w::CCASensor) = pressure(w.calibr)
caltemp(w::CCASensor) = reftemp(w.calibr)
 
function velocity(w::CCASensor, E;
                  T=caltemp(w), P=pressure(w),
                  fluid=fluid(w), I=current(w))
    g = gain(w)
    Ew = E / g # Voltage passing through the resistor
    Rw = Ew / I # Resistance of the sensor
    return velocity(w.calibr, E; T=T, P=P, fluid=fluid, I=I)
        
end

function velocity!(U::AbstractArray, w::CCASensor, E::AbstractArray;
                   T=caltemp(w), P=pressure(w),
                   fluid=fluid(w), I=current(w))
    U .= velocity.(w, E; T=T, P=P, fluid=fluid, I=I)
    return U
end

velocity(w::CCASensor, E::AbstractArray;
         T=caltemp(w), P=pressure(w),
         fluid=fluid(w), I=current(w)) =
             velocity!(similar(E), w, E, T=T, P=P, fluid=fluid, I=I)

                    

(w::CCASensor)(args...; kw...) = velocity(w, args...; kw...)
