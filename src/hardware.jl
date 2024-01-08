
abstract type AbstractHardware end

hardmodel(h::AbstractHardware) = h.model
hardtag(h::AbstractHardware) = h.tag

struct HWCable{U,V,W} <: AbstractHardware
    "Model/serial/etc of the cable"
    model::String
    "Tag"
    tag::String
    "Resistance of the cable"
    R::U
    "Impedance of the cable"
    impedance::V
    "Length of the cable"
    L::W
end

HWCable(R,impedance; model="", tag="") = HWCable(code, tag, R, impedance)


impedance(c::HWCable) = c.impedance

resistance(c::HWCable) = resistance(c.R)
resistance(c::HWCable, T) = resistance(c.R, T)





struct HWSupport{U} <: AbstractHardware
    "Support model"
    model::String
    "Support tag (storage and control)"
    tag::String
    "Nominal resistance of the support"
    R::U
end

"""
    `HWSupport(;model="55H21", tag="", R=0.44, D=4.0, L=235, Lc=765)`

Stores information about Probe support. Noty strictly necessary but 
could be useful for future reference

"""
HWSupport(R ; model="", tag="") = HWSupport(model, tag, R)

resistance(c::HWSupport) = resistance(c.R)
resistance(c::HWSupport, T) = resistance(c.R, T)




