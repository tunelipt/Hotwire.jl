module Hotwire

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export CTASensor, Wire, resistor, optemperature, overtemp, overheat_ratio, gain
export CCASensor, current
export tempcorr
include("resistor.jl")
include("wire.jl")
include("tempcorr.jl")

include("probes.jl")
end # module
