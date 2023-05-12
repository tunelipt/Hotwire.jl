# Temperature/Pressure/etc correction

abstract type AbstractAnemometerCorrection end

abstract type AbstractCTACorrection <: AbstractAnemometerCorrection
abstract type AbstractCCACorrection <: AbstractAnemometerCorrection


struct CTACorrection <: AbstractCCACorrection

end

