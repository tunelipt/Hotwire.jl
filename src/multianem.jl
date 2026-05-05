# multianem.jl
# Systems with the same circuit but several sensors (assumed to be the same).

struct Multianem{Anem<:AbstractThermalAnemometer,Setup,Device}
    sensor::Vector{Anem}
    setup::Setup
    dev::Device
end


function tomltomultiname(toml)
    # The offchans tells us what channels should not be used
    if haskey(toml, "offchans")
        offchans = Int.(toml["offchans"])
    else
        offchans = Int[]
    end
    

    # First read the resistance elements
    # If it is a single element, then the same resistance is used for every
    # sensor.
    if !haskey(toml, "resistance")
        error("A multianemometer should have at least 1 resistance element")
    end
    tres = toml["resistance"]
    if tres isa AbstractVector
        # We will assume that the same resistance type for every resistance element
        r0 = tomltoresistance(tres[1])
        R = typeof(r0)[r0]
        for t in tres[begin+1:end]
            push!(R, tomltoresistance(t))
        end
        nwires = length(tres)
    else
        # We should specify the number of wires in nwires
        nwires = toml["nwires"]
        R = fill(tomltoresistance(tres), nwires)
    end

    # Now read the bridge information
    tbridge = toml["bridge"]
    if tbridge isa AbstractVector
        # Each bridge is specified independently
        # They should have the same type
        # The number should bw the same as the number of resisance elements:
        length(tbridge) != nwires && error("The number of bridges should be the same as sensors in a `multianem` anemometer")

        bridges = [tomltobridge(tb) for tb in tbridge]
    else
        bridges = fill(tomltobridge(tbridge), nwires)
    end
    # The bridges should have the same mode:
    modes = unique(b.mode for b in bridges)
    length(modes) > 1  && error("All bridges should have the same mode, 'cta', or 'cca'")
    mode = modes[1]
    
    # Now we should read the calibration curves
    # There might not be any calibration and in this case
    # there will no correction model 
    if !haskey(toml, "calibration")
        calibr = [NoCalibr(R[1], [0.0], [0.0], 1.0, 0.0, 0.0, x->x; fluid=AIR)
                  for i in 1:nwires]
        nocal = true
    else
        # Read calibration of every sensor
        cal = [tomltocalibration(c) for c in toml["calibration"]]

        # Read the Correction scheme
        # Usually it will be only one but it could be more.
        # In this case it should have the same number of wires
        tcorr = toml["correction"]
        if tcorr isa AbstractVector
            length(tcorr) != nwires && error("The number of correction sections should be the same as the number of wires if an array is wanted")
            correct = [tomltocorrection(t) for t in tcorr]
        else
            correct = fill(tomltocorrection(tcorr), nwires)
        end
        nocal = false
        
    end

    # Read the data acquisition info
    # Simplest case, an array of [[daq]] or single device and an array of channels
    if !haskey(toml, "daq")
        # No data acquisition
        daq = nothing
    else
        tdaq = toml["daq"]
        if tdaq isa AbstractVector
            length(tdaq) != nwires && error("daq is an array and should have the same length as `nwires`")
            daq = [tomltodaq(t) for t in tdaq]
        else
            daq1 = tomltodaq(tdaq)
            length(daq1.chan) != nwires && error("Incompatible number of sensors and daq channels")
            daq = [(dev=daq.dev, chan=c) for c in daq1.chan]
        end
    end

    # Let's build the data structures
    idx = trues(nwires)
    idx[offchans] .= false

    
    
end


