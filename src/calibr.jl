export velocity, velocity!, HultmarkSmits2010


abstract type AbstractAnemCalibr end

caltemp(cal::AbstractAnemCalibr) = cal.T
calpress(cal::AbstractAnemCalibr) = cal.P
calwtemp(cal::AbstractAnemCalibr) = cal.Tw
calwres(cal::AbstractAnemCalibr) = cal.Rw
calfluid(cal::AbstractAnemCalibr) = cal.fluid


function makecaltable(E::AbstractVector{X}, U::AbstractVector{X}, T, P) where {X}
    N = length(E)
    @assert size(E,1) == size(U,1)

    caltab = zeros(T,N,4)

    caltab[:,1] .= E
    caltab[:,2] .= U
    caltab[:,3] .= T
    caltab[:,4] .= P

    Tm = mean(T)
    Pm = mean(P)
    
    return Tm, Pm, caltab
end


struct TempCalibr{X,RT<:AbstractResistor,Fluid,Fit} <: AbstractAnemCalibr
    "Resistance element"
    R::RT
    "Calibration operating resistance"
    Rw::X
    "Calibration operating resistance temperature"
    Tw::X
    "Fluid calibration temperature"
    T::X
    "Fluid calibration pressure"
    P::X
    "Calibration fluid"
    fluid::Fluid
    "Calibration data"
    data::Matrix{X}
    "Calibration curve fit"
    fit::Fit
end

Base.broadcastable(calibr::AbstractAnemCalibr) = Ref(calibr)


function TempCalibr(R::RT,
                    Ec::AbstractVector,
                    Uc::AbstractVector,
                    Tc, Pc, Rwc, makefitfun;
                    Rw=nothing, T=nothing, P=nothing,
                    fluid=AIR) where {RT<:AbstractResistor}

    Tm, Pm, caltab = makecaltable(Ec, Uc, Tc, Pc)

    # Reference values
    if isnothing(T)
        T = Tm
    end
    if isnothing(P)
        P = Pm
    end
    if isnothing(Rw)
        Rw = mean(Rwc)
    end

    Twc = temperature.(R, Rwc)
    
    Tw = temperature(R, Rw)

    # Correct voltages to conditions Rw, T, P
    E = sqrt.( Rw./Rwc .* (Tw - T) ./ (Twc .- Tc)  )
    
    fit = makefitfun(E, U)

    
    return TempCalibr(R, Rw, Tw, T, P, fluid, caltab, fit)
    
end

function hwcorrect(cal::TempCalibr; T=caltemp(cal), P=calpress(P),
                   fluid=calfluid(cal),Rw=calwres(cal))
    
    Tw = temperature(cal.R, Rw)
    return sqrt( cal.Rw/Rw * (cal.Tw - cal.T) / (Tw - T) )
end

    
function velocity(cal::TempCalibr, E; T=caltemp(cal), P=calpress(P),
                  fluid=calfluid(cal),Rw=calwres(cal))
    fcorr = hwcorrect(cal; T=T, P=P, fluid=fluid, Rw=Rw)
    Ec = fcorr * E
    return fit(E)
end

function velocity!(U::AbstractArray, cal::TempCalibr, E::AbstractArray;
                   T=caltemp(cal), P=calpress(P),
                   fluid=calfluid(cal),Rw=calwres(cal))
    @assert size(U) == size(E)
    fcorr = hwcorrect(cal; T=T, P=P, fluid=flud, Rw=Rw)
    for i in eachindex(E)
        U[i] = cal.fit(E[i]  * fcorr)
    end
    return U
        
end

struct HWCalibr{X<:AbstractFloat,RT<:AbstractResistor,
                         Fluid,Fit} <: AbstractAnemCalibr
    "Resistance element"
    R::RT
    "Calibration operating resistance"
    Rw::X
    "Calibration operating resistance temperature"
    Tw::X
    "Fluid calibration temperature"
    T::X
    "Fluid calibration pressure"
    P::X
    "Calibration fluid"
    fluid::Fluid
    "Calibration data"
    data::Matrix{X}
    "Calibration curve fit"
    fit::Fit
    "Prandtl number exponent in Nusselt"
    n::X
    "Mean surface temperature factor"
    theta::X
    "Calibration kinematic viscosity"
    nu::X
end


function HWCalibr(R::RT,
                           Ec::AbstractVector,
                           Uc::AbstractVector,
                           Tc, Pc, Rwc, makefitfun; n=0.3, theta=0.5,
                           Rw=nothing, T=nothing, P=nothing,
                           fluid=AIR) where {RT<:AbstractResistor}

    Tm, Pm, caltab = makecaltable(Ec, Uc, Tc, Pc)
    
    # Reference values
    if isnothing(T)
        T = Tm
    end
    if isnothing(P)
        P = Pm
    end
    if isnothing(Rw)
        Rw = mean(Rwc)
    end
    
    Twc = temperature.(R, Rwc)
    
    Tw = temperature(R, Rw)
    ΔTc = Twc .- Tc
    ΔT  = Tw - T
    
    # Film temperature
    Tfc = Tc .+ theta .* ΔTc  # Individual calibration points
    Tf  = T   + theta  * ΔT  # Calibration reference conditions
    
    phic   = heatcond.(fluid, Tfc,  Pc) .* prandtl.(fluid, Tfc, Pc) .^ n

    # Reynolds number
    nuc = kinvisc.(fluid, Tfc, Pc)
    nu  = kinvisc(fluid, Tf, Pc)
    
    Rec = Uc ./ nuc

    # E²/ϕΔT
    Ex = Ec .* Ec ./ (phic .* (Twc .- Tc))
   

    fit = makefitfun(Ex, Rec)

   
    
    return HWCalibr(R, Rw, Tw, T, P, fluid, caltab, fit, n, nu)
    
end



function velocity(cal::HWCalibr, E; T=caltemp(cal), P=calpress(P),
                  fluid=calfluid(cal),Rw=calwres(cal))

    Tw = temperature(cal.R, Rw)

    ΔT = Tw - T
    Tf = T + cal.theta * ΔT
    
    k = heatcond(fluid, Tf, P)
    Pr = prandtl(fluid, Tf, P)
    phi = k*Pr^cal.n

    nu = kinvisc(fluid, Tf, P)
    
    return nu * cal.fit(E*E/(phi * Rw * ΔT))
end


function velocity!(U::AbstractArray, cal::HWCalibr, E::AbstractArray;
                   T=caltemp(cal), P=calpress(P),
                   fluid=calfluid(cal),Rw=calwres(cal))

    @assert size(E) == size(U)
    
    Tw = temperature(cal.R, Rw)
    ΔT = Tw - T
    Tf = T + cal.theta * ΔT
    k = heatcond(fluid, Tf, P)
    Pr = prandtl(fluid, Tf, P)
    den = k*Pr^cal.n * Rw * (Tw - T)
    
    nu = kinvisc(fluid, Tf, P)

    for (e,i) in enumerate(E)
        U[i] = nu * cal.fit(e*e/den)
    end
    return U
end



