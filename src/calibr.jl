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


(cal::CalibrCurve)(E) = cal.fit(E)



