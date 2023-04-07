module Hotwire

export Thermistor, Resistor, AbstractResistor, temperature, resistance, reftemp
export CTASensor, Wire, resistor, optemperature, overtemp, overheat_ratio, gain
export AbstractThermalAnemometer
export CCASensor, current
export tempcorr
export calibr_curve, CalibrCurve, ExtrapolateFit
export HWKingFit, BlendFit

export AbstractProbe, AbstractProbe1d, AbstractProbe2d, AbstractProbe3d
export HWCable, HWSupport, HWBridge
export Probe1d, Probe2d, Probe3d
export velocity, velocity!
export dircalibr
include("resistor.jl")
include("wire.jl")
include("tempcorr.jl")
include("fit.jl")
include("calibr.jl")
include("probes.jl")
include("probe2d.jl")
include("probe3d.jl")
include("hardware.jl")

end # module
