using Statistics


abstract type AbstractCalibr end
abstract type AbstractCalibr1d end

struct CalibrCurve{Correct,Fluid,Fit,U} <: AbstractCalibr1d
    "Anemometer calibration curve"
    fit::Fit
    "Calibration fluid"
    fluid::Fluid
    "Operating conditions of the anemometer"
    temp::TempCorrect{U}
    "Correction algorithm"
    corr::Correct
end

Base.broadcastable(cal::CalibrCurve) = Ref(cal)

reftemp(cal::CalibrCurve) = reftemp(cal.temp)
resistance(cal::CalibrCurve) = resistance(cal.temp)
temperature(cal::CalibrCurve) = temperature(cal.temp)
fluid(cal::CalibrCurve) = cal.fluid


function correction(cal::CalibrCurve; E, T=reftemp(cal), P=101325.0,
                    Tw=temperature(cal), Rw=resiscance(cal),
                    fluid=fluid(cal))
    
    
end



