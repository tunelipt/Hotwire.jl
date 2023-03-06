struct HWSupport
    "Support model"
    model::String
    "Support tag (storage and control)"
    tag::String
    "Nominal resistance of the support"
    R::Float64
    "Diameter of the support"
    D::Float64
    "Length of the support"
    L::Float64
    "Length of support cables"
    Lc::Float64
end

"""
    `HWSupport(;model="55H21", tag="", R=0.44, D=4.0, L=235, Lc=765)`

Stores information about Probe support. Noty strictly necessary but 
could be useful for future reference

"""
HWSupport(;model="", tag="", R=0.0, D=0.0, L=0.0, Lc=0.0) = HWSupport(model, tag, R, D, L, Lc)

struct HWCable
    "Cable model/type"
    model::String
    "Cable tag (storage and control)"
    tag::String
    "Nominal resistance of the cable"
    R::Float64
    "Length of the cable"
    L::Float64
    "Connector types (left and right)"
    conn::Tuple{String,String}
end

"""
    `HWCable(model="", tag="", R=0.2, L=4.0, ("BNC","BNC"))`

Stores information about the cable used. Not strictly necessary but
could be useful.

"""
HWCable(;model="A1863 - 50Ω impedance", tag="", R=0.2, L=4.0, conn=("BNC", "BNC")) = HWCable(model, tag, R, L, conn)

struct HWBridge
    "Anemometer system model"
    model::String
    "Anemometer tag"
    tag::String
    "Specific bridge id"
    id::Int
    "Output offset"
    offset::Float64
    "Output gain"
    gain::Float64
    "Output low pass filter"
    lpfilt::Float64
    "Bridge ration"
    bridgeratio::String
    "Bridge amplifier gain"
    ampgain::Float64
    "Bridge filter"
    filter::Float64
    "Cable compensation"
    cablecomp::Int
    "Cooling interval"
    coolinterv::Float64
end

"""
    `HWBridge(model="", tag="", id=1, offset=0, gain=1, lpfilt=1e3,
              bridgeration="20:1", ampgain=8, filter=3, 
              cablecomp=0, coolinterv=2)`

Information about hotwire hardware system.

"""
HWBridge(;model="Streamline", tag="", id=1, offset=0.0, gain=1.0,
         lpfilt=1e3, bridgeratio="20:1", ampgain=8, filter=3,
         cablecomp=0, coolinterv=2.0) =
             HWBridge(model, tag, id, offset, gain, lpfilt,
                      bridgeratio, ampgain, filter, cablecomp,
                      coolinterv)
