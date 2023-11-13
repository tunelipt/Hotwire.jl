using Statistics


abstract type AbstractCalibr end
abstract type AbstractCalibr1d end

struct CalibrCurve{Fit,Tab} <: AbstractCalibr1d
    "Anemometer calibration curve"
    fit::Fit
    "Anemometer calibration data"
    data::Tab
end
Base.broadcastable(cal::CalibrCurve) = Ref(cal)


(cal::CalibrCurve)(E) = cal.fit(E)



