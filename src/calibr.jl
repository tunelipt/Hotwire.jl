using Statistics


abstract type AbstractCalibr end
abstract type AbstractCalibr1d end

"""
`CalibrCurve(fit, fluid, temp_corr, corr)`

Creates a `CalibrCurve` object that handles the calibration curve of a
thermal anemometer.

## Arguments

 * `fit` an object that, given a voltage, returns the velocity in calibration conditions
 * `fluid` The fluid that was used as calibration and is the default fluid when using
 * `temp` An object, usually [`TempCorrect`](@ref) that stores calibration conditions of
the anemometer: reference temperature, operating temperature and operating resistance. See below.
 * `corr` A correction model. See below

When using constant current anemometers, during calibration as the fluid velocity
increases, the sensor temperature will decrease. But for `CalibrCurve`, these changing
values are corrected for constant resistance and temperature.

Different models to correct are possible. This package implements the following models
 * [`TOnlyCorrect`](@ref) only operating temperatures and resistance are considered
 * [`WireCorrect`](@ref) uses convective heat transfer concepts to correc sensors whose surface temperature corresponds to the temperature of the resistive element of the sensor.
 * [`GlassbeadCorrect`](@ref) Implements a model where the resistive element has an insulating shell and wire leads.
 * [`mf58correct`](@ref) A specific implementation of the [`GlassbeadCorrect`](@ref) model

"""
struct CalibrCurve{Correct,Fluid,Fit,TCorr,U} <: AbstractCalibr1d
    "Anemometer calibration curve"
    fit::Fit
    "Calibration fluid"
    fluid::Fluid
    "Operating conditions of the anemometer"
    temp::TCorr
    "Correction algorithm"
    corr::Correct
    "Calibration pressure in Pa"
    P::U
end

Base.broadcastable(cal::CalibrCurve) = Ref(cal)

reftemp(cal::CalibrCurve) = reftemp(cal.temp)
resistance(cal::CalibrCurve) = resistance(cal.temp)
temperature(cal::CalibrCurve) = temperature(cal.temp)
fluid(cal::CalibrCurve) = cal.fluid
pressure(cal::CalibrCurve) = cal.P

correctmodel(cal::CalibrCurve) = cal.corr



"""
`correction(cal, E, tc, mc)`

Corrects anemometer output of an anemoeter to different operating conditions

## Arguments
 * `cal` [`CalibrCurve`](@ref) object with calibration of the sensor
 * `E` voltage accross resistance element of the sensor
 * `tc` [`TempCorrect`](@ref) with operating conditions
 * `mc` [`AbstractAnemCorrect`](@ref) correction model in operating conditions



"""
function correction(cal::CalibrCurve{Correct}, E, 
                    tc::TempCorrect, mc::Correct) where {Correct}
    
    anemcorrect(E, cal.temp, cal.corr, tc, mc)
    
end

function correction(cal::CalibrCurve{Correct}, E, tc::TempCorrect,
                    fluid=fluid(cal), P=pressure(cal)) where {Correct}
    
    anemcorrect(E, cal.temp, cal.corr, tc, mc)
    
end


function correction(cal::CalibrCurve{Correct}, E, tc, fluid, P)
    mc_cal = 


(cal::CalibrCurve)(E::Real) = cal.fit(E)



