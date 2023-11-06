
abstract type AbstractThermalAnemometer end
abstract type AbstractConstCurrAnemometer <: AbstractThermalAnemometer end


"""
    CCASensor(R, a)

A structure to manage constant current anemometer sensors (CCA)

 * `R` Resistance element
 * `I` current in Ampere


"""
struct CCASensor{T,U,V,ResType<:AbstractResistor} <: AbstractConstCurrAnemometer
    "Element whose resistance changes with temperature"
    R::ResType
    "Operating current in Ampere"
    I::T
    "Output gain"
    gain::U
    "Reference temperature"
    T₀::V
end
CCASensor(R::ResType, I, T0=reftemp(R); gain=1) where {ResType<:AbstractResistor} =
    CCASensor(R, I, gain, T0)
Base.broadcastable(sensor::CCASensor) = Ref(sensor)

"Current in A flowing through the CCA"
current(w::CCASensor) = w.I

reftemp(w::CCASensor) = w.T₀
