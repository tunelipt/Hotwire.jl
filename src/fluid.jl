
import Polynomials: ImmutablePolynomial

const IPoly = ImmutablePolynomial


const Ru = 8314.46261815324 # J/(K.kmol)

const M_Air   = 28.96518
const M_N₂    = 28.01340
const M_C₃H₈  = 44.09562
const M_CO₂   = 44.00980 
const M_He    =  4.00205
const M_O₂    = 31.99880
const M_H₂    =  2.01588

"""
`ConstPropFluid(rho, mu, k, cp)`

Creates a fluid with constant properties. The user should provide
the values of density, dynamic viscosity, thermal conductivity
and specific heat, all in SI units.
"""
struct ConstPropFluid{T,U,V,W}
    "Density of fluid in kg/m³"
    ρ::T
    "Dynamic Viscosity of fluid in Pa⋅s"
    μ::U
    "Thermal conductivity of fluid in W/m⋅K"
    k::V
    "Prandtl Number of fluid"
    cp::W
end
Base.broadcastable(fluid::ConstPropFluid) = Ref(fluid)

ConstPropFluid(fluid, T, P) = ConstPropFluid(density(fluid, T, P),
                                             viscosity(fluid, T, P),
                                             heatcond(fluid, T, P),
                                             specheat(fluid, T, P))

"Computes the heat conductivity of a fluid in W/mK"
heatcond(x::ConstPropFluid, T, P=101325.0) = x.k

"Computes the dynamic viscosity of a fluid in Pa⋅s"
viscosity(x::ConstPropFluid, T, P=101325.0) = x.μ

"Computes the density of a fluid in kg/m³"
density(x::ConstPropFluid, T, P=101325.0) = x.ρ

"COmputes the specific heat of a fluid in J/Kg⋅K"
specheat(x::ConstPropFluid, T, P=101325.0) = x.cp





"Computes the Prandtl number of a fluid"
prandtl(x, T, P=101325.0) = specheat(x,T,P) * viscosity(x,T,P) / heatcond(x,T,P)

"Computes the kinematic viscosity of a fluid in m²/s"
kinvisc(x, T, P=101325.0) = viscosity(x,T,P) / density(x,T,P)


"""
`ConstPropFluidProp(ρ, μ, k, Pr)
`ConstPropFluidProp(;rho=ρ, mu=μ, k=k, Pr=Pr)`

Models a fluid with constant properties.

"""
ConstPropFluid(;rho=1.1, mu=1.8e-5, k=26e-3, cp=1005.0) =
    ConstPropFluid(rho, mu, k, cp)


"""
    `IdealGas(M, c_cp, c_mu, c_k)`

Creates an `IdealGas` object that is used to compute thermodynamic and transport
properties of ideal gases. The main parameter is the mean molecular mass of the gas
in kg/kmol.

Functions for calculating the specific heat at constant pressure, dynamic
viscosity and 
"""
struct IdealGas{U,FV,FW,FX} 
    M::U
    cpfun::FV
    mufun::FW
    kfun::FX
end
Base.broadcastable(fluid::IdealGas) = Ref(fluid)



specheat(x::IdealGas, T, P=101325.0) = x.cpfun(T) 
viscosity(x::IdealGas, T, P=101325.0) = x.mufun(T)
heatcond(x::IdealGas, T, P=101325.0) = x.kfun(T)
density(x::IdealGas, T, P=101325.0) = P * x.M / (Ru*T)








const c_cₚ_C₃H₈ = Ru/M_C₃H₈ .* (4.21093013E+00, 1.70886504E-03, 7.06530164E-05,
                               -9.20060565E-08, 3.64618453E-11)
const c_μ_C₃H₈  = (-2.62714e-7, 2.87582e-8, 8.02292e-13, -1.00381e-14, 3.67025e-18)
const c_k_C₃H₈  = (-0.00190848, 2.62236e-5, 1.37295e-7, 8.1589e-12, -2.96982e-15)


const c_cₚ_Air = Ru/M_Air .* (3.56839620E+00, -6.78729429E-04, 1.55371476E-06,
                             -3.29937060E-12, -4.66395387E-13)
const c_μ_Air = (-9.87746e-8, 7.85222e-8, -6.67915e-11, 4.13539e-14)
const c_k_Air = (5.88994e-5, 9.54173e-5, -2.39694e-8, -1.12551e-11, 4.85253e-15)


const c_cₚ_N₂ = Ru/M_N₂ .* (3.53100528E+00, -1.23660988E-04, -5.02999433E-07,
                          2.43530612E-09, -1.40881235E-12)
const c_μ_N₂ = (-5.26692e-7, 8.28046e-8, -9.36255e-11, 8.18452e-14, -3.1254e-17,
                -1.01669e-17)
const c_k_N₂ = (-0.0006241, 0.000109906, -9.13581e-8, 7.26363e-11, -2.73799e-14)

const c_cₚ_CO₂ = Ru/M_CO₂ .* (0.23568130E+01, 0.89841299E-02,-0.71220632E-05,
                             0.24573008E-08,-0.14288548E-12)

const c_μ_CO₂ = (7.19059e-6, 5.43706e-10, 1.22949e-10, -1.50777e-13, 5.9986e-17)
const c_k_CO₂ = (0.00522826, -1.43376e-5, 2.41587e-7, -2.49604e-10, 8.42276e-14)


const c_cₚ_He  = (Ru / M_He * 2.50000000E+00, )
const c_μ_He = (4.09517e-6, 5.55975e-8, 2.90184e-12, -5.28468e-14, 3.6016e-17)
const c_k_He = (0.0301597, 0.000452246, -1.56078e-7, 2.31013e-11)

const c_cₚ_O₂ = Ru / M_O₂ .* (3.78245636E+00, -2.99673416E-03,  9.84730201E-06, 
                             -9.68129509E-09, 3.24372837E-12)
const c_μ_O₂ = (-9.66187e-7, 9.56013e-8, -9.94813e-11, 8.15837e-14, -2.95515e-17)
const c_k_O₂ = (-0.0018068, 0.000123112, -1.45337e-7, 1.83738e-10, -9.52034e-14)


const c_cₚ_H₂  = Ru / M_H₂ .* (2.34433112E+00, 7.98052075E-03,-1.94781510E-05,
                              2.01572094E-08, -7.37611761E-12)

const c_k_H₂ = (-0.0155079, 0.000945333, -1.33176e-6, 1.43404e-9, -5.97372e-13)
const c_μ_H₂ = (1.04885e-6, 3.57237e-8, -4.36667e-11, 4.7382e-14, -2.12645e-17)

       
const AIR = IdealGas(M_Air, IPoly(c_cₚ_Air,:T), IPoly(c_μ_Air,:T), IPoly(c_k_Air,:T))
const NITROGEN = IdealGas(M_N₂, IPoly(c_cₚ_N₂,:T), IPoly(c_μ_N₂,:T), IPoly(c_k_N₂,:T))
const OXYGEN = IdealGas(M_O₂,  IPoly(c_cₚ_O₂,:T), IPoly(c_μ_O₂,:T), IPoly(c_k_O₂,:T))
const HELIUM = IdealGas(M_He, IPoly(c_cₚ_He,:T), IPoly(c_μ_He,:T), IPoly(c_k_He,:T))
const C3H8 = IdealGas(M_C₃H₈, IPoly(c_cₚ_C₃H₈,:T), IPoly(c_μ_C₃H₈,:T),
                      IPoly(c_k_C₃H₈,:T))

const CO2  = IdealGas(M_CO₂, IPoly(c_cₚ_CO₂,:T), IPoly(c_μ_CO₂,:T), IPoly(c_k_CO₂,:T))

const HYDROGEN  = IdealGas(M_H₂, IPoly(c_cₚ_H₂,:T), IPoly(c_μ_H₂,:T), IPoly(c_k_H₂,:T))


    
