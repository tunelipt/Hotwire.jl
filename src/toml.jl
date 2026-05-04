# Configuring anemometers from TOML

function toml_hotwire(toml)
    t = toml["type"] 
    if t != "hotwire"
        error("Type expected: `hotwire` got `$t`")
    end
    
    if !haskey(toml, "nwires") || toml["nwires"] == 1
        # We will assume 1 wire
        return toml_hw1d(toml)
    else
        Nw = toml["nwires"]
    end

        
end


function toml_hw1d(toml)
    
    # First we need to read the resistance
    if !haskey(toml, "resistance")
        error("A thermal anemometer probe should have a resistance element!")
    else
        R = toml_resistance(toml["resistance"])
    end

    # Now we read the bridge (hardware) configuration
    if !haskey(toml, "bridge")
        error("An anemometer should specify a bridge (hardware configuration)")
    else
        bridge = toml_bridge(toml["bridge"])
    end

    # This is the data acquisition
    if haskey(toml, "daq")
        daq = toml_daq(toml["daq"])
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
            cal = toml_calibration(toml["calibration"])
        end
    end

    # Now the correction model
    if !isnothing(cal)
        if !haskey(toml, "correction")
            error("With a calibration, there should be a correction model!")
        end
        correction = toml_correction(toml["correction"])
        
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

function toml_bridge(toml)

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



function toml_resistance(toml)
    ks = keys(toml)
    if "type" ∉ ks
        # If not provided we will assume linear resistor
        t = "linear"
    else
        t = toml["type"]
    end

    if t == "linear"
        R0 = toml["R"]
        a  = toml["a"]
        T0 = toml["T"]
        R = Resistor(R=R0, a=a, T=T0)
    elseif t == "ntc" # Negtative temperature coefficient thermistor
        R0 = toml["R"]
        B  = toml["B"]
        T0 = toml["T"]
        R = Thermistor(R=R0, B=B, T=T0)
    else
        error("Unknown resistance element type `$t`")
    end

    return R
        
end


function toml_fluid(toml)
    if toml isa AbstractString
        s = string(toml)
        if !haskey(standard_ideal_fluid_table, s)
            error("Fluid $s is not currently available")
        else
            return standard_ideal_fluid_table[s]
        end
    else
        if !haskey(toml, "type")
            # We should specify somehow the type of fluid model here
            error("I don't know how to create the fluid. You should specify the type")
        else
            t = toml["type"]
            if t == "constant"
                # We should specify the following properties:
                # rho, mu, k, cp.
                rho = toml["rho"]  # Density
                mu = toml["mu"]  # Dynamic viscosity
                k  = toml["k"]   # thermal conductivity function
                cp = toml["cp"]  # Specific heat at constant pressure
                return ConstPropFluid(rho, mu, k, cp)
            elseif t == "coolprop_constant"
                fl = toml["composition"]
                T = toml["T"]
                P = toml["P"]
                rho = PropsSI("D", "P", P, "T", T, fl)
                mu = PropsSI("V", "P", P, "T", T, fl)
                k = PropsSI("L", "P", P, "T", T, fl)
                cp = PropsSI("C", "P", P, "T", T, fl)
                return ConstPropFluid(rho, mu, k, cp)
            elseif t == "coolprop"
                # For now only pure fluids should be used
                # I will try to improve on that later. CoolProp is really cool...
                return CPFluid(toml["composition"])
            elseif t == "humidair"
                # We will also use CoolProp
                # But the Psychrometrics stuff.
                T = 293.15
                P = 101325.0
                hvar = "R"  # Relative humidity
                hval  = 0.5  # 50%
                humvars = ["W", "B", "D", "R", "Y"] 
                for (kk,vv) in toml
                    
                    if kk == "T"
                        T = vv
                    elseif kk == "P"
                        P = vv
                    elseif kk ∈ humvars
                        hvar = kk
                        hval = vv
                    end
                end
                return HumidAir(T, P, hvar, hval)
            else
                error("Unknown fluid type $t")
            end
        end
        
    end
end


function toml_calibration(toml)

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
        fluid = toml_fluid(toml["fluid"])
    end

    return (U=U, E=E, T=T, P=P, fluid=fluid)
    
end


function toml_fit(toml)
    model = toml["model"]

    if model == "KingPoly"
        n = toml["n"]
        a = toml["a"]
        fun = (x,y) -> KingPoly(x, y, a, n)
    elseif model == "Polynomial"
        model == "Polynomial"
        n = toml["n"]
        fun =  (x,y) -> Polynomials.fit(x, y, n)
    else
        error("Fit model is $model. Unknown, choose `KingPoly` or `Polynomial`")
    end

    return fun
end

function toml_correction(toml)
    model = toml["model"]

    if model == "HWCalibr"
        model1 = HWCalibr
        params = (n = toml["n"], theta = toml["theta"])
        fit = toml_fit(toml["fit"])
    elseif model == "TempCalibr" || model == "NoCorrection"
        model1 = TempCalibr
        params = nothing
        fit = toml_fit(toml["fit"])
    elseif model == "NoCorrection"
        model1 = NoCorrection
        params = ()
        fit = toml_fit(toml["fit"])
    else
        error("Model $model  not implemented!")
    end

    return (model=model, fun=model1, params=params, fit=fit)
end

function toml_daq(toml)
    device = toml["device"]
    channel = toml["channel"]

    return(dev=device, chan=channel)
end
