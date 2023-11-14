module Hotwire

import Polynomials
import Polynomials: Polynomial

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export reftemp, refresist
export CTASensor, Wire, resistor, overheatratio, overtemp, gain
export AbstractThermalAnemometer
export CCASensor, current
export tempcorr
export calibr_curve, CalibrCurve
export KingLaw, KingPoly, Polynomial
export make_king_fitfun, make_kingpoly_fitfun, make_poly_fitfun

export AbstractProbe, AbstractProbe1d, AbstractProbe2d, AbstractProbe3d
export HWCable, HWSupport, HWBridge
export Probe1d, Probe2d, Probe3d
export velocity, velocity!
export dircalibr


abstract type AbstractThermalAnemometer end
abstract type AbstractCTA <: AbstractThermalAnemometer end
abstract type AbstractCCA <: AbstractThermalAnemometer end

"Resistor element of a thermal anemometer sensor"
resistor(w::AbstractThermalAnemometer) = w.R

"Gain of anemometer output"
gain(w::AbstractThermalAnemometer) = w.gain


include("resistor.jl")
include("correct.jl")
include("ctasensor.jl")
#include("wire.jl")
include("fluid.jl")
#include("tempcorr.jl")
#include("fit.jl")
#include("calibr.jl")
#include("probes.jl")
#include("probe2d.jl")
#include("probe3d.jl")
#include("hardware.jl")

end # module
