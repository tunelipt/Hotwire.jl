abstract type AbstractAnemCalibr end

caltemp(cal::AbstractAnemCalibr) = cal.T
calpress(cal::AbstractAnemCalibr) = cal.P
calwtemp(cal::AbstractAnemCalibr) = cal.Tw
calwres(cal::AbstractAnemCalibr) = cal.Rw
calfluid(cal::AbstractAnemCalibr) = cal.fluid


function velocity(cal::AbstractAnemCalibr, E; T=caltemp(cal), P=calpress(P),
                  fluid=calfluid(cal),Tw=calwtemp(cal), Rw=calwres(cal))
    corr = correct(cal, T=T, P=P, fluid=fluid, Tw=Tw, Rw=Rw)
    return velocity(cal, E, corr)
end

function velocity!(U::AbstractArray, cal::AbstractAnemCalibr, E::AbstractArray, corr)
    for i in eachindex(E)
        U[i] = velocity(cal, E[i], corr)
    end
    return U
end

function velocity!(U::AbstractArray, cal::AbstractAnemCalibr, E::AbstractArray;
                   T=caltemp(cal), P=calpress(P),
                   fluid=calfluid(cal),Tw=calwtemp(cal), Rw=calwres(cal))
    corr = correct(cal, T=T, P=P, fluid=fluid, Tw=Tw, Rw=Rw)
    for i in eachindex(E)
        U[i] = velocity(cal, E[i], corr)
    end
    return U
end

velocity(cal::AbstractAnemCalibr, E::AbstractArray, corr) = velocity!(similar(E), cal,
                                                                      E, corr)


velocity(cal::AbstractAnemCalibr, E::AbstractArray;
         T=caltemp(cal), P=calpress(P),
         fluid=calfluid(cal),Tw=calwtemp(cal), Rw=calwres(cal)) =
             velocity!(similar(E), cal, E, correct(cal; T=T, P=P, Tw=Tw, Rw=Rw,
                                                   fluid=fluid))



struct TempCalibr{T,Fluid,Fit} <: AbstractAnemCalibr
    "Calibration operating resistance"
    Rw::T
    "Calibration operating resistance temperature"
    Tw::T
    "Fluid calibration temperature"
    T::T
    "Fluid calibration pressure"
    P::T
    "Calibration fluid"
    fluid::Fluid
    "Calibration curve fit"
    fit::Fit
    "Prandtl number exponent in Nusselt"
    n::T
    "k.Prⁿ at calibration conditions"
    phi::T
    "Calibration kinematic viscosity"
    nu::T
end



struct HultSmits2010{T<:AbstractFloat,Fluid,Fit} <: AbstractAnemCalibr
    "Calibration operating resistance"
    Rw::T
    "Calibration operating resistance temperature"
    Tw::T
    "Fluid calibration temperature"
    T::T
    "Fluid calibration pressure"
    P::T
    "Calibration fluid"
    fluid::Fluid
    "Calibration curve fit"
    fit::Fit
    "Prandtl number exponent in Nusselt"
    n::T
    "k.Prⁿ at calibration conditions"
    phi::T
    "Calibration kinematic viscosity"
    nu::T
end

#=
function calcorrect(::Type{HultSmits2010}, E, T, P, Tw, Rw, fluid,
                    Tc, Pc, Twc, Rwc; n=0.3)
    Tf   = (Tw + T) / 2
    ΔTw  = Tw - T

    Pr = prandtl(fluid, Tf, P)
    k  = heatcond(fluid, Tf, P)
    nu = kinvisc(fluid, Tf, P)
    phi = k * Pr^n

    # Reference calibration conditions
    Tfc   = (Twc + Tc) / 2
    ΔTwc  = Twc - Tc

    Prc = prandtl(fluid, Tfc, Pc)
    kc  = heatcond(fluid, Tfc, Pc)
    nuc = kinvisc(fluid, Tfc, Pc)
    phic = kc * Prc^n
=#
    
   
    

function HultSmits(R::RT, E::AbstractVector, U::AbstractVector, T, P,
                   Rw, makefitfun;
                   Rwc=nothing, Tc=nothing, Pc=nothing,
                   fluid=AIR, n=0.3) where {RT<:AbstractResistor}

    Tw = temperature.(R, Rw)

    # Default calibration parameters will be the mean parameters during calibration
    isnothing(Rwc) && Rwc = mean(Tw)  #
    Twc = temperature(R, Rwc)
    isnothing(Tc) && Tc = mean(T) #
    isnothing(Pc) && Pc = mean(P)  
    
    # Properties at reference calibration conditions (Rwc, Twc, Pc, Tc)
    Tfc = (Tc + Twc) / 2
    Prc = prandtl(fluid, Tfc, Pc)
    kc  = heatcond(fluid, Tfc, Pc)
    nuc = kinvisc(fluid, Tfc, Pc)
    
    phic = kc * Prc^n
    ΔTwc = Tw - T

    # Conditions during calibration
    Tf = 
    ef = E .^ 2 ./ (phi*Rw*ΔTw)
    uf = U ./ nu
    fit = makefitfun(ef, uf)
    return HultSmits2010(Rw, Tw, T, P, fluid, fit, n, phi, nu)
end


function correct(cal::HultSmits2010; T=caltemp(cal), P=calpress(P),
                 fluid=calfluid(cal),Tw=calwtemp(cal), Rw=calwres(cal))
    Tf = (T + Tw) / 2

    k = heatcond(fluid, Tf, P)
    Pr = prandtl(fluid, Tf, P)
    phi = k*Pr^cal.n
    nu = kinvisc(fluid, Tf, P)
    
    ef = (cal.phi/phi) * (cal.Rw / Rw) * (cal.Tw - cal.T) / (Tw - T)
    uf = nu 
    return (ef, uf)
end


function velocity(cal::HultSmits2010, E, corr)
    ef = corr[1] * (E*E / (cal.phi * cal.Rw * (cal.Tw - cal.T)))
    return corr[2] * cal.fit(ef)
end



struct HWCalibr{U,RT<:AbstractResistor,Fluid,Fit} <: AbstractAnemCalibr
    "Resistance element"
    R::RT
    "Calibration operating resistance"
    Rw::U
    "Calibration operating resistance temperature"
    Tw::U
    "Fluid calibration temperature"
    T::U
    "Fluid calibration pressure"
    P::U
    "Calibration fluid"
    fluid::Fluid
    "Calibration curve fit"
    fit::Fit
    "Prandtl number exponent in Nusselt"
    n::U
    "Reynolds number exponent in Nusselt"
    m::U
    "Nondimensional film temperature"
    theta::U
    "Calibration (k.Prⁿ)^(1/m) / ν"
    phi::U
end


function HWCalibr(E::AbstractVector{X}, U::AbstractVector{X}, makefitfun;
                  Tw, Rw, T, P, fluid=AIR, m=0.45, n=0.3, theta=0.5) where {X}

    ΔTw = Tw - T
    Tf = T + theta * ΔTw # Film temperature

    fe = 
    
end

