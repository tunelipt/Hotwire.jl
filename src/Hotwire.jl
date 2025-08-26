module Hotwire

import Polynomials
import Polynomials: Polynomial
import Statistics: mean

export Thermistor, Resistor, AbstractResistor, temperature, resistance
export reftemp, refresist

export ConstPropFluid, heatcond, prandtl, viscosity, kinvisc, density, specheat
export AIR, HELIUM, NITROGEN, C3H8, OXYGEN, IdealGas, CO2, HYDROGEN

export TempCalibr, HWCalibr 

export AbstractAnemCalibr, calpress, caltemp, calwtemp, calwres, calfluid
export AbstractAnemCorrect, TempCorrect, WireCorrect, InsulatedCorrect
export mf58correct, mf52correct
export correct, tempcorrect
export CorrFactor, corrfactor, voltage, sensorvolt, outsignal, linsignal

export CalibrCurve, correctmodel, correction, pressure

export CTASensor, Wire, resistor, overheatratio, overtemp, gain,caltemp
export AbstractThermalAnemometer
export CCASensor, current
export KingLaw, KingPoly, Polynomial, PowerPoly
export make_king_fitfun, make_kingpoly_fitfun, make_poly_fitfun

export sensor
export AbstractProbe, AbstractProbe1d, AbstractProbe2d, AbstractProbe3d
export HWCable, HWSupport, HWBridge, impedance
export Probe1d, Probe2d, Probe3d
export velocity, velocity!, velf, velf!
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
include("calibr.jl")
include("ctasensor.jl")
#include("ccasensor.jl")
include("fit.jl")
#include("probes.jl")
#include("probe2d.jl")
#include("probe3d.jl")
#include("hardware.jl")
#include("dantec.jl")

end # module
