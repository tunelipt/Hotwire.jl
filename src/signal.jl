
# This file implements a relationship between resistor voltage
# and signal output.

struct LinearTransform{U,V}
    gain::U
    offset::V
end

linsignal(g,o) = LinearTransform(g,o)

"""
`sensorvolt(s, E)`

Returns the voltage accross the sensor element taking into consideration
the output voltage that has a gain and an offset.
"""
sensorvolt(w,E) = E/w.gain + w.offset

"""
`outsignal(w,E)`

Returns the output voltage by applying sensor gain and offset
from voltage accross the sensor.
"""
outsignal(w, E) = (E - w.offset)*w.gain


applycorr(w, E, fc::CorrFactor) = outsignal(sensorvolt(w,E)*fc.f)

# Valid only for `LinearTransform` but useful since this is the
# most probable case...
gain(w) = w.gain
offset(w) = w.offset
    
    

