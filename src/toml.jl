# Configuring anemometers from TOML

function tomltohotwire(toml)
    t = toml["type"] 
    if t != "hotwire"
        error("Type expected: `hotwire` got `$t`")
    end
    
    if !haskey(toml, "nwires") || toml["nwires"] == 1
        # We will assume 1 wire
        return tomltohw1d(toml)
    else
        Nw = toml["nwires"]
    end

        
end


function tomltohw1d(toml)
    
    # First we need to read the resistance
    if !haskey(toml, "resistance")
        error("A thermal anemometer probe should have a resistance element!")
    else
        R = tomltoresistance(toml["resistance"])
    end

    # Now we read the bridge (hardware) configuration
    if !haskey(toml, "bridge")
        error("An anemometer should specify a bridge (hardware configuration)")
    else
        bridge = tomltobridge(toml["bridge"])
    end

    # This is the data acquisition
    if haskey(toml, "daq")
        daq = tomltodaq(toml["daq"])
    else
        daq = nothing
    end

    if !haskey(toml, "calibration")
        error("If you just want the voltages, no calibration stuff, set `calibration = false`")
    else
        if !(toml isa AbstractDict)
            if toml["calibration"] != false
                error("If there is no calibration table, `calibration` should be `false`")
            else
                cal = nothing
            end
        else
            cal = tomltocalibration(toml["calibration"])
        end
    end

    # Now the correction model
    if !isnothing(cal)
        if !haskey(toml, "correction")
            error("With a calibration, there should be a correction model!")
        end
        correction = tomltocorrection(toml["correction"])
        
    else
        correction = nothing
    end

    # Let's build the probe
    
    if bridge.mode == :cta
        sensor = build_cta(R, bridge, cal, correction)
    elseif bridge.mode == :cca
        sensor = build_cca(R, bridge, cal, correction)
    end
    
    return Probe1d(sensor, toml, daq)
        
end

function build_cta(R, bridge, cal, correction)
    Rw = bridge.params.R

    if isnothing(cal)
        calibr = NoCalibr(293.15, 101325.0, AIR)
    else
        calibr = correction.fun(R, cal.U, cal.E, Rw,
                                cal.T, cal.P, correction.fit;
                                fluid=cal.fluid, correction.params...)
    end
    return CTASensor(R, Rw, calibr)
        
end

function build_cca(R, bridge, cal, correction)
    I = bridge.params.I
    if isnothing(cal)
        calibr = NoCalibr(293.15, 101325.0, AIR)
    else
        Rw = cal.E ./ I
        calibr = correction.fun(R, cal.U, cal.E, Rw,
                                cal.T, cal.P, correction.fit;
                                fluid=cal.fluid, correction.params...)
    end
    return CCASensor(R, I, calibr)
        
end

function tomltobridge(toml)

    # Check the mode
    if !haskey(toml, "mode")
        m = :cta
    else
        m = Symbol(lowercase(toml["mode"]))
        if m != :cta && m != :cca
            error("Hardware mode $m should be `cta` or `cca`")
        end
    end

                     
    if m == :cta
        # We should read the resistance
        if !haskey(toml, "R")
            error("A `cta` probe should specify the operating resistance")
        else
            params = (R=toml["R"],)
        end
                     
    else
        if !haskey(toml, "I")
            error("A `cca` probe should specify the operating current")
        else
            params = (I=toml["I"],)
        end
    end

    # Now we will read the signal processing stuff
    if haskey(toml, "gain")
        gain = toml["gain"]
    else
        gain = 1.0
    end
    
    if haskey(toml, "offset")
        offset = toml["offset"]
    else
        offset = 0.0
    end

    return (mode=m, params=params, gain=gain, offset=offset)
end







function tomltocalibration(toml)

    # A calibration has the fields E, U, T, P and fluid
    if !haskey(toml, "E")
        error("Probe calibration should have field `E` with voltage")
    else
        E = Float64.(toml["E"])
    end
    if !haskey(toml, "U")
        error("Probe calibration should have field `U` with velocity")
    else
        U = Float64.(toml["U"])
    end
    Ne = length(E)
    Nu = length(U)
    @assert Ne == Nu "U ($Nu elements) should have the same length as E ($Ne elements)"
    
    if !haskey(toml, "T")
        error("Probe calibration should have field `T` with temperature")
    else
        T1 = toml["T"]
        if T1 isa AbstractVector
            @assert length(T1)==Ne "`T` should have a single value or the same number of elements of `U`, $Nu elements"
            T = Float64.(T1)
        else
            T = fill(T1, Nu)
        end
    end

    if !haskey(toml, "P")
        error("Probe calibration should have field `T` with pressure")
    else
        P1 = toml["P"]
        if P1 isa AbstractVector
            @assert length(T1)==Ne "`P` should have a single value or the same number of elements of `U`, $Nu elements"
            P = Float64.(P1)
        else
            P = fill(P1, Nu)
        end
    end

    # Now we should read the fluid. If not provided we will assume AIR
    if !haskey(toml, "fluid")
        fluid = AIR
    else
        fluid = tomltofluid(toml["fluid"])
    end

    return (U=U, E=E, T=T, P=P, fluid=fluid)
    
end



function tomltocorrection(toml)
    model = toml["model"]

    if model == "HWCalibr"
        model1 = HWCalibr
        params = (n = toml["n"], theta = toml["theta"])
        fit = tomltofit(toml["fit"])
    elseif model == "TempCalibr" || model == "NoCorrection"
        model1 = TempCalibr
        params = nothing
        fit = tomltofit(toml["fit"])
    elseif model == "NoCorrection"
        model1 = NoCorrection
        params = ()
        fit = tomltofit(toml["fit"])
    else
        error("Model $model  not implemented!")
    end

    return (model=model, fun=model1, params=params, fit=fit)
end

function tomltodaq(toml)
    device = toml["device"]
    channel = toml["channel"]

    return(dev=device, chan=channel)
end
