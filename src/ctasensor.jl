

struct CTASensor{T,U,RT,Calibr,Correct} <: AbstractCTA
    R::RT
    Rw::T
    gain::U
    cal::Calibr
    corr::Correct
end

Base.broadcastable(sensor::CTASensor) = Ref(sensor)


sensor(w::AbstractCTA) = w.R
resistance(w::AbstractCTA) = w.Rw
temperature(w::AbstractCTA) = temperature(sensor(w), resistance(w))
overheatratio(w::AbstractCTA, T) = resistance(w) / resistance(sensor(w),T) - 1
overheatratio(w::AbstractCTA) = resistance(w) / resistance(sensor(w)) - 1
overtemp(w::AbstractCTA,T) = temperature(sensor(w),resistance(w)) - T
overtemp(w::AbstractCTA) = temperature(sensor(w),resistance(w)) - temperature(sensor(w))



