module Hotwire

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export CTASensor, resistor, optemperature, overtemp, overheat_ratio, gain
export CCASensor, current
export tempcorr
include("resistor.jl")
include("wire.jl")
include("tempcorr.jl")

end # module
