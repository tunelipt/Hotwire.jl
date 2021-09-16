module Hotwire

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export CTASensor, Wire, resistor, optemperature, overtemp, overheat_ratio, gain
export AbstractThermalAnemometer
export CCASensor, current
export tempcorr
export CalibrCurve, ExtrapolateFit

include("resistor.jl")
include("wire.jl")
include("tempcorr.jl")
include("probes.jl")
include("calibr.jl")
end # module
