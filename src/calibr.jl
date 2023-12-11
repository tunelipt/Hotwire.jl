using Statistics


abstract type AbstractCalibr end
abstract type AbstractCalibr1d end

struct CalibrCurve{Correct,Fluid,Fit,TCorr} <: AbstractCalibr1d
    "Anemometer calibration curve"
    fit::Fit
    "Calibration fluid"
    fluid::Fluid
    "Operating conditions of the anemometer"
    temp::TCorr
    "Correction algorithm"
    corr::Correct
end

Base.broadcastable(cal::CalibrCurve) = Ref(cal)

reftemp(cal::CalibrCurve) = reftemp(cal.temp)
resistance(cal::CalibrCurve) = resistance(cal.temp)
temperature(cal::CalibrCurve) = temperature(cal.temp)
fluid(cal::CalibrCurve) = cal.fluid



function correction(E, cal::CalibrCurve{Correct},
                    tc::TempCorr, mc::Correct) where {Correct}
    
    anemcorrect(E, cal.temp, cal.corr, tc, mc)
    
end


function correctmodel(cal::CalibrCurve{Correct},
                      tc::TempCorr, fluid, P=101325.0) where {Correct}

    Correct(cal.corr, tc, fluid, P)
end

function velocity(E, cal::CalibrCurve, tc::TempCorr, cm::AbstractAnemCorrect)
    Ec = anemcorrect(E, cal.temp, cal.corr, tc, cm)
    Uc = cal.fit(Ec)
    return kinvisc(cm) / kinvisc(cal.corr)
end

(cal::CalibrCurve)(E::Real) = cal.fit(E)

function (cal::CalibrCurve{AC})(E::Real, tc::TempCorr,
                                mc::AC) where {AC<:AbstractAnemCorrect}
end


