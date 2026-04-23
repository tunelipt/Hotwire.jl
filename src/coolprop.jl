import CoolProp: PropsSI, HAPropsSI

export CPFluid, HumidAir

struct CPFluid
    fluid::String
end
Base.broadcastable(fluid::CPFluid) = Ref(fluid)

heatcond(fl::CPFluid, T, P=101325.0) = PropsSI("L", "P", P, "T", T, fl.fluid)
viscosity(fl::CPFluid, T, P=101325.0) = PropsSI("V", "P", P, "T", T, fl.fluid)
density(fl::CPFluid, T, P=101325.0) = PropsSI("D", "P", P, "T", T, fl.fluid)
specheat(fl::CPFluid, T, P=101325.0) = PropsSI("C", "P", P, "T", T, fl.fluid)
prandtl(fl::CPFluid, T, P=101325.0) = PropsSI("PRANDTL", "P", P, "T", T, fl.fluid)


# Humid Air, we are going to use transport properties for dry air since they
# don't appear to be implemented on coolprop

struct HumidAir
    w::Float64
end

function HumidAir(;T=293.15, P=101325.0, hum...)
    @assert length(hum) == 1 error("$hum should have a single argument specifying the humidity")
    h = string(first(keys(hum)))
    hval = first(values(hum))
    
    w = HAPropsSI("W", "T", T, "P", P, h, hval)
    HumidAir(w)
end

specheat(fl::HumidAir, T, P=101325.0) = HAPropsSI("cp_ha", "P", P, "T", T, "W", fl.w)
heatcond(fl::HumidAir, T, P=101325.0) = PropsSI("L", "P", P, "T", T, "Air")
density(fl::HumidAir, T, P=101325.0) = 1/HAPropsSI("Vha", "P", P, "T", T, "W", fl.w)
viscosity(fl::HumidAir, T, P=101325.0) = PropsSI("V", "P", P, "T", T, "Air")
prandtl(fl::HumidAir, T, P=101325.0) = PropsSI("PRANDTL", "P", P, "T", T, "Air")
