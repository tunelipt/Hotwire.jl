# DANTEC probes and hardware

struct StreamlineBridge <: AbstractHardware
    "Anemometer system model"
    model::String
    "Anemometer tag"
    tag::String
    "Overheat ratio"
    overheat::Float64
    "Decade resistance"
    Rdecade::Float64
    "Probe resistance (Ω)"
    Rprobe::Float64
    "Cable resistance (Ω)"
    Rcable::Float64
    "Sensor resistance (Ω)"
    R0::Float64
    "Reference temperature (K)"
    Tr::Float64
    "Ambient Temperature (K)"
    Ta::Float64
    "Bridge top voltage (V)"
    brdgtop::Float64
    "Low pass filter (kHz)"
    lpfilt::Float64
    "High pass filter (kHz)"
    hpfilt::Float64
    "Signal offset (V)"
    offset::Float64
    "Signal gain"
    gain::Float64
    "Bridge resistance ratio"
    bridgeratio::Float64
    "Bridge filter"
    filter::Int16
    "Bridge gain"
    bgain::Int16
    "Cable compensation"
    cablecomp::Int16
    "Specific bridge type"
    type::Int16
end

"""
`StreamlineBridge(config; model="", tag="", T0=nothing, br=1/20)`

Information about hotwire hardware system.

"""
function StreamlineBridge(config ;model="Streamline", tag="", T0=nothing, br=1/20)
    btype = round(Int16,config[1])
    overheat = confgig[2]
    Rdecade = config[3]
    Rprobe = config[4]
    Rcable = config[5]
    R0 = config[6]
    if !isnothing(T0)
        Tr = T0
        Ta = T0
    else
        Tr = config[7]
        Ta = config[8]
    end
    brdgtop = config[9]
    lpfilt = config[10]
    hpfilt = config[11]
    offset = config[12]
    gain = config[13]
    bridgeratio = br
    filter = round(Int16,config[14])
    bgain =  round(Int16,config[15])
    cablecomp =  round(Int16,config[16])

    StreamlineBridge(model, tag, overheat, Rdecade, Rprobe, Rcable, R0, Tr, Ta,
                     brdgtop, lpfilt, hpfilt, offset, gain, bridgeratio,
                     filter, bgain, cablecomp, btype)
end
    

function dantec1d(hconfig, calibr, fitfun; model="", tag="", alpha=0.4e-2, fluid=AIR,
                  T0=nothing, br=1/20, probe=nothing, support=nothing, cable=nothing,
                  n=1/3)
    
    bridge = StreamlineBridge(hconfig, model=model, tag=tag, T0=T0, br=br)

    # Let's calculate the operating temperature
    Rdecade   = bridge.Rdecade
    Rprobe    = bridge.Rprobe
    R0        = bridge.R0

    ΔR =  (Rdecade * br  - Rprobe)
    Rw = R0 + ΔR

    T0 = bridge.Tr + 273.15  # Resistor reference temperature in K

    R = Resistor(R=R0, a=alpha, T=T0)
    Tw = temperature(R, Rw)
    
    Tcal = calibr[:,3] .+ 273.15  # Calibration temperature (K)
    Pcal = calibr[:,4] .* 1000    # Calibration pressure (Pa)
    E    = calibr[:,2]
    Uc   = calibr[:,1]
    
    # Mean temperature and pressure - they will be used
    Tm = mean(Tcal)
    Pm = mean(Pcal)

    g = bridge.gain
    o = bridge.offset
    
    corr = WireCorrect(Tm, Pm, fluid, Rw, Rw, n)
    
    Ei = [(e/g + o) for e in E]
    fc = [correct(e, corr, Tm, Pm, fluid, Rw, Tw) for e in Ei]
    Ec = [(f.E - o)*g for f in fc]

    fit = fitfun(Ec, Uc)

    sensor = CTASensor(R, Rw, Tw, g, o, corr, fit)
    
    return Probe1d(sensor, (probe, support, cable, bridge))
end
