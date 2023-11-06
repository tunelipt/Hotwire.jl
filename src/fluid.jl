export AbstractFluidProp, ConstPropFluid, IdealGas, Air


const Ru = 8314.46261815324 # J/(K.kmol)

abstract type AbstractFluidProp end

struct ConstPropFluid{T,U,V,W} <: AbstractFluidProp
    "Density of fluid in kg/m³"
    ρ::T
    "Dynamic Viscosity of fluid in Pa⋅s"
    μ::U
    "Thermal conductivity of fluid in W/m⋅K"
    k::V
    "Prandtl Number of fluid"
    Pr::W
end

"""
`ConstPropFluidProp(ρ, μ, k, Pr)
`ConstPropFluidProp(;rho=ρ, mu=μ, k=k, Pr=Pr)`

Models a fluid with constant properties.

"""
ConstPropFluid(;rho=1.1, mu=1.8e-5, k=26e-3, Pr=0.7) =
    ConstPropFluid(rho, mu, k, Pr)

fluidprops(x::ConstPropFluid, T, P) = x.ρ, x.μ, x.k, x.Pr


struct IdealGas{U,V} <: AbstractFluidProp
    R::U
    propfun::V
end

function IdealGas(;mw=28.96518, mu=1.813e-5, k=25.87e-3, Pr=0.707)
    R = Ru/mw
    IdealGas(Rg, (T,P)->(P/(R*T), mu, k, Pr))
end


fluidprops(x::IdealGas, T, P) = x.propfun(x.R, T, P)

function airvisc(T,P)
    S = 113
    μ₁ = 18.13e-6
    T₁ = 293.15
    μ₁ * (T/T₁)^1.5 * (T₁ + S) / (T + S)
end

const Mair = 28.96518


const Ckair = [-1.433633e-4, 1.0184e-4, -4.8574e-8, 1.5207e-11]
const Cpair = [3.56839620E+00, -6.78729429E-04, 1.55371476E-06, 
               -3.29937060E-12, -4.66395387E-13]

aircond(T,P) = evalpoly(T, Ckair) + 0

    
airspecheat(T,P) = evalpoly(T, Cpair) * Ru/Mair



function airprops(R, T, P)
    ρ = P / (R*T)
    μ = airvisc(T,P)
    k = aircond(T,P)
    cₚ = airspecheat(T,P) 
    Pr = cₚ * μ / k
    return ρ,μ,k,Pr
end



#=
Burcat data for air

132259-10-0
AIR calculated from ingredients %N2=78.084 %O2=20.9476 %Ar=0.9365 %CO2=0.0319
This format is not capable of automatic formula calculation for this species!!!
See New NASA Polynomials REF=McBride & Gordon NASA RP-1271 1992 Sructure
for automatic formula calculation should be:
N 1.56O 0.42AR 0.01C 0.00 Max Lst Sq Error Cp @ 2500 K 0.19%
AIR L 9/95 WARNING! 0.G 200.000 6000.000 B 28.96518 1
3.08792717E+00 1.24597184E-03-4.23718945E-07 6.74774789E-11-3.97076972E-15 2
-9.95262755E+02 5.95960930E+00 3.56839620E+00-6.78729429E-04 1.55371476E-06 3
-3.29937060E-12-4.66395387E-13-1.06234659E+03 3.71582965E+00-1.50965000E+01 4
=#

       
const Air = IdealGas(Ru/Mair, airprops)







    
