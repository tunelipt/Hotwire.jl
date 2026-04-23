# Calibration table

struct CalibrTable{X<:AbstractFloat,Fluid}
    "Calibration velocity"
    U::Vector{X}
    "Calibration temperature"
    T::Vector{X}
    "Calibration pressure"
    P::Vector{X}
    "Calibration output for each wire"
    E::Matrix{X}
    "Calibration fluid"
    fluid::Fluid
    "Wire order"
    wires::Vector{Int}
end
