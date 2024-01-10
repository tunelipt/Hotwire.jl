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

    signal = linsignal(bridge.gain, bridge.offset)
    
    corr = WireCorrect(Tm, Pm, fluid, Rw, Rw, n)
    
    Ei = [(e/g + o) for e in E]
    fc = [correct(e, corr, Tm, Pm, fluid, Rw, Tw) for e in Ei]
    Ec = [(f.E - o)*g for f in fc]

    fit = fitfun(Ec, Uc)

    sensor = CTASensor(R, Rw, Tw, signal, corr, fit)
    
    return Probe1d(sensor, (hconfig, calibr, probe, support, cable, bridge))
end


function dantec2d(hconfig, calibr, fitfun, k²;
                  model="", tag="", alpha=0.4e-2, fluid=AIR,
                  T0=nothing, br=1/20, probe=nothing, support=nothing, cable=nothing,
                  n=1/3)
    
    bridge = [StreamlineBridge(hconfig[i,:], model=model, tag=tag, T0=T0, br=br)
              for i in 1:2]
    
    T0 = [b.Tr for b in bridge] .+ 273.15
    Rdecade   = [b.Rdecade for b in bridge]
    Rprobe    = [b.Rprobe for b in bridge]
    R0        = [b.R0 for b in bridge]
    

    ΔR =  (Rdecade * br  - Rprobe)
    Rw = R0 + ΔR


    R = [Resistor(R=R0[i], a=alpha, T=T0[i]) for i in eachindex(bridge)]
    Tw = temperature.(R, Rw)

    Tcal = calibr[:,4] .+ 273.15  # Calibration temperature (K)
    Pcal = calibr[:,5] .* 1000    # Calibration pressure (Pa)
    E1    = calibr[:,2]
    E2    = calibr[:,3]
    Uc    = calibr[:,1]

    
    Tm = mean(Tcal)
    Pm = mean(Pcal)

    g = [b.gain for b in bridge]
    o = [b.offset for b in bridge]
    signal1 = linsignal(g[1], o[1])
    signal2 = linsignal(g[2], o[2])
    
    corr = WireCorrect(Tm, Pm, fluid, Rw, Rw, n)

    Ei1 = [(e/g[1] + o[1]) for e in E1]
    Ei2 = [(e/g[2] + o[2]) for e in E2]

    fc1 = [correct(e, corr, Tcal[i], Pcal[i], fluid, Rw, Tw)
           for (i,e) in enumerate(Ei1)]
    fc2 = [correct(e, corr, Tcal[i], Pcal[i], fluid, Rw, Tw)
           for (i,e) in enumerate(Ei2)]

    Ec1 = [outsignal(signal1, e*f.f) for (e,f) in zip(fc1,Ei1)]
    Ec2 = [outsignal(signal2, e*f.f) for (e,f) in zip(fc2,Ei2)]

    fit1 = fitfun(Ec1, Uc)
    fit2 = fitfun(Ec2, Uc)

    
    sensor1 = CTASensor(R[1], Rw[1], Tw[1], signal1, corr, fit1)
    sensor2 = CTASensor(R[2], Rw[2], Tw[2], signal2, corr, fit2)

    return Probe2d((sensor1, sensor2), (k²[1], k²[2]),
                   (hconfig, calibr, probe, support, cable, bridge))
    
end


function dantec3d(hconfig, calibr, fitfun, k², h²;
                  ih=[3,1,2], idircal=1, cosang=cosϕ_dantec,
                  model="", tag="", alpha=0.4e-2, fluid=AIR,T0=nothing,
                  br=1/20, probe=nothing, support=nothing, cable=nothing,
                  n=1/3)
    
    bridge = [StreamlineBridge(hconfig[i,:],
                               model=model, tag=tag, T0=T0, br=br)
              for i in 1:3]
    
    T0 = [b.Tr for b in bridge] .+ 273.15
    Rdecade   = [b.Rdecade for b in bridge]
    Rprobe    = [b.Rprobe for b in bridge]
    R0        = [b.R0 for b in bridge]
    

    ΔR =  (Rdecade * br  - Rprobe)
    Rw = R0 + ΔR


    R = [Resistor(R=R0[i], a=alpha, T=T0[i]) for i in eachindex(bridge)]
    Tw = temperature.(R, Rw)

    Tcal = calibr[:,5] .+ 273.15  # Calibration temperature (K)
    Pcal = calibr[:,6] .* 1000    # Calibration pressure (Pa)
    E1    = calibr[:,2]
    E2    = calibr[:,3]
    E3    = calibr[:,4]
    Uc    = calibr[:,1]
    
    Tm = mean(Tcal)
    Pm = mean(Pcal)

    g = [b.gain for b in bridge]
    o = [b.offset for b in bridge]
    signal1 = linsignal.(g[1], o[1])
    signal2 = linsignal.(g[2], o[2])
    signal3 = linsignal.(g[3], o[3])
    
    corr = WireCorrect(Tm, Pm, fluid, Rw, Rw, n)

    Ei1 = [sensorvolt(signal[1], e) for e in E1]
    Ei2 = [sensorvolt(signal[2], e) for e in E2]
    Ei3 = [sensorvolt(signal[3], e) for e in E3]

    # 
    fc1 = [correct(e, corr, Tcal[i], Pcal[i], fluid, Rw, Tw)
           for (i,e) in enumerate(Ei1)]
    fc2 = [correct(e, corr, Tcal[i], Pcal[i], fluid, Rw, Tw)
           for (i,e) in enumerate(Ei2)]
    fc3 = [correct(e, corr, Tcal[i], Pcal[i], fluid, Rw, Tw)
           for (i,e) in enumerate(Ei3)]
    
    Ec1 = [outsignal(signal1, e*f.f) for (e,f) in zip(fc1,Ei1)]
    Ec2 = [outsignal(signal1, e*f.f) for (e,f) in zip(fc2,Ei2)]
    Ec3 = [outsignal(signal1, e*f.f) for (e,f) in zip(fc3,Ei3)]

    fit1 = fitfun(Ec1, Uc)
    fit2 = fitfun(Ec2, Uc)
    fit3 = fitfun(Ec3, Uc)

    
    sensor1 = CTASensor(R[1], Rw[1], Tw[1], signal1, corr, fit1)
    sensor2 = CTASensor(R[2], Rw[2], Tw[2], signal2, corr, fit2)
    sensor3 = CTASensor(R[3], Rw[3], Tw[3], signal3, corr, fit3)

    return Probe3d((sensor1, sensor2, signal3), (k²[1], k²[2], k²[3]),
                   (h²[1], h²[2], h²[3]), 
                   (hconfig, calibr, probe, support, cable, bridge);
                   ih=ih, idircal=idircal, cosang=cosang)
    
    
end
