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


function tomlconfig_fluid(toml)
    if toml isa AbstractString
        s = string(toml)
        if s ∉ keys(standard_ideal_fluid_table)
            error("Fluid $s is not currently available")
        else
            return standard_ideal_fluid_table[s]
        end
    else
        if "type" ∉ keys(toml)
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
