module Hotwire

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export CTASensor, resistor, optemperature, overtemp, overheat_ratio
include("resistor.jl")
include("wire.jl")


end # module
