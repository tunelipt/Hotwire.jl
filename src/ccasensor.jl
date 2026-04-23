

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
struct CCASensor{U,RT,Calibr} <: AbstractCCA
    "Temperature dependent resistor"
    R::RT
    "Operating electrical current"
    I::U
    "Output gain"
    gain::U
    "Calibration and correction model"
    calibr::Calibr
end
Base.broadcastable(sensor::CCASensor) = Ref(sensor)

resistor(w::AbstractCCA) = w.R

"Operating current of the CTA"
current(w::AbstractCCA) = w.Rw

gain(w::AbstractCCA) = w.gain
reftemp(w::CCASensor) = reftemp(w.R)
refresist(w::CCASensor) = refresist(w.R)


calibration(w::AbstractCCA) = w.calibr


 
function velocity(w::CCASensor, E;
                  T=caltemp(w), P=calpress(w),
                  fluid=calfluid(w), I=current(w))
    g = gain(w)
    Ew = E / g # Voltage passing through the resistor
    Rw = Ew / I # Resistance of the sensor
    return velocity(calibration(w), resistor(w), E, Rw,  T, P, fluid)
        
end

function velocity!(U::AbstractArray, w::CCASensor, E::AbstractArray;
                   T=caltemp(w), P=pressure(w),
                   fluid=fluid(w), I=current(w))

    R = resistor(w)
    g = gain(w)
    cal = calibration(w)
    map!(e->velocity(cal, R, e, resistance(R,e/(g*I)), T, P, fluid), U, E)
end

function velocity(w::CCASensor, E::AbstractArray;
                  T=caltemp(w), P=pressure(w),
                  fluid=fluid(w), I=current(w))
    velocity!(similar(E), w, E; T=T, P=P, fluid=fluid, I=I)
end

(w::CCASensor)(args...; kw...) = velocity(w, args...; kw...)
