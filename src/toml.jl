# Configuring anemometers from TOML

function hwconfig(toml)
    ks = keys(toml)
    t = toml["type"] 
    if t != "hotwire"
        error("Type expected: `hotwire` got `$t`")
    end

    if "nwires" ∉ ks
        # We will assume 1 wire
        return tomlconfig1d(toml)
    else
        Nw = toml["nwires"]
    end

    
        
        
end

function tomlconfig1d(toml)
    
    # the mode of operation CTA or CCA
    ks = keys(toml)
    if "cta" ∈ ks
        # We need the operating resistance
        if "R" ∉ toml["cta"]
            error("A CTA should have the operating resistance specified as field `R`")
        else
            Rw = toml["cta"]["R"]
        end
    end
    
    # We should get the resistance element
    if "resistance" ∉ ks
        error("A resistance element should be present in `resistance` field")
    else
        R = tomlconfig_resistance(toml["resistance"])
    end
    
end

function tomlconfig_resistance(toml)
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
