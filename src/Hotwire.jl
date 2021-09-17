module Hotwire

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export CTASensor, Wire, resistor, optemperature, overtemp, overheat_ratio, gain
export AbstractThermalAnemometer
export CCASensor, current
export tempcorr
export calibr_curve, CalibrCurve, ExtrapolateFit

export AbstractProbe, AbstractProbe1d, AbstractProbe2d, AbstractProbe3d
export HWCable, HWSupport, HWBridge
export Probe1d, Probe2d, Probe3d


include("resistor.jl")
include("wire.jl")
include("tempcorr.jl")
include("probes.jl")
include("calibr.jl")
end # module
