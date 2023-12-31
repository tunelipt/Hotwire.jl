module Hotwire

import Polynomials
import Polynomials: Polynomial

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export reftemp, refresist

export ConstPropFluid, heatcond, prandtl, viscosity, kinvisc, density, specheat
export AIR, HELIUM, NITROGEN, C3H8, OXYGEN, IdealGas

export AbstractAnemCorrect, TempCorrect, WireCorrect, GlassbeadCorrect
export mf58correct, correct, tempcorrect
    

export CalibrCurve, correctmodel, correction, pressure

export CTASensor, Wire, resistor, overheatratio, overtemp, gain
export AbstractThermalAnemometer
export CCASensor, current
export KingLaw, KingPoly, Polynomial
export make_king_fitfun, make_kingpoly_fitfun, make_poly_fitfun

export sensor
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
include("fluid.jl")
include("correct.jl")
include("ctasensor.jl")
include("ccasensor.jl")
include("fit.jl")
#include("probes.jl")
#include("probe2d.jl")
#include("probe3d.jl")
#include("hardware.jl")

end # module
