
import Statistics: mean
"""
`example_calibr()`

Returns a typical data from a calibration curve from o constant temperature hotwire.


"""
function example_calibr(fluid=AIR)
    A = 1.45
    B = 0.6
    n = 0.43
    king = (A=A, B=B, n=n)

    Tc = 20.0 + 273.15
    Pc = 101325.0

    Rdec  =  120.0
    Rtot  =    3.850 # Total probe resistance
    Rcs   =    0.650 # cable+support resistance
    Rl    =    0.50  # lead resistance
    R0    = Rtot - Rcs - Rl  # Sensor cold resistance
    T0    =  293.15  # Reference temperature
    br = 1/20  # Bridge ratio
    ΔR = Rdec * br - Rtot
    Rw = R0 + ΔR
    a = ΔR / R0
    
    R = Resistor(R=R0, a=0.45e-2, T=T0)
    
    U = [0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0, 15.0, 18.0]

    E = sqrt.(A .+ B .* U.^n)
    Tw = temperature(R, Rw)
    return (E=E, U=U, T=Tc, P=Pc, Rw=Rw, Tw=Tw, fluid=fluid, R=R, king=king)
end


let
    X = example_calibr(AIR) #ConstPropFluid)
    king = X.king
    # Let's create an exact equation for the calibration curve
    fit = E -> ( (E*E - king.A) / king.B )^(1/king.n)
    ifit = U -> king.A + king.B * U ^ king.n
    
    Tw = temperature(X.R, X.Rw)
    cal = TempCalibr(X.R, X.U, X.E, X.Rw, X.T, X.P, makekingfitfun(a=king.n, N=1),
                      fluid=AIR)
    @test isapprox(X.T, cal.T)
    @test isapprox(X.P, cal.P)
    

    # Test variations in Operating resistance
    # We will create a new calibration where Rw is changed by a fator of f.
    f = 1.1
    α = X.R.a
    calx = TempCalibr(X.R, X.U, X.E, f*X.Rw, X.T, X.P, makekingfitfun(a=king.n, N=1),
                      fluid=AIR)
        
    @test all(isapprox.(velocity.(cal, X.R,  X.E, X.Rw,   X.T, X.P, X.fluid),
                        velocity.(calx, X.R, X.E, X.Rw*f, X.T, X.P, X.fluid)))

    

    # Changes in pressure: we should not have any!
    @test all(isapprox.(velocity.(cal, X.R, X.E, X.Rw, X.T, X.P, X.fluid),
                        velocity.(cal, X.R, X.E, X.Rw, X.T, X.P*0.9, X.fluid)))

    # Check changes in temperature
    Ta = cal.T + 10.0
    E1 = X.E .* sqrt.( (X.Tw - Ta) / (X.Tw - cal.T) )
    U1 = velocity.(cal, X.R, E1, X.Rw, Ta, X.P, X.fluid)
    U2 = velocity.(cal, X.R, X.E, X.Rw, X.T, X.P, X.fluid)
    @test all(isapprox.(U1, U2))

    
    
end



# Now we will test HWCalibr - using Constant properties!
let
    n = 0.3
    theta = 0.5
    
    fluid  = ConstPropFluid(AIR, 293.15, 101325.0)

    X = example_calibr(fluid) #ConstPropFluid)
    king = X.king
    # Let's create an exact equation for the calibration curve
    fit = E -> ( (E*E - king.A) / king.B )^(1/king.n)
    ifit = U -> king.A + king.B * U ^ king.n

    
    # Correct only temperature. Since properties are constant,
    # it should be identical to TempCalibr
    calt = TempCalibr(X.R, X.U, X.E, X.Rw, X.T, X.P, makekingfitfun(a=king.n, N=1),
                      fluid=fluid)
    cal = HWCalibr(X.R, X.U, X.E, X.Rw, X.T, X.P, makekingfitfun(a=king.n, N=1),
                   fluid=fluid, n=n, theta=theta)
    @test all(isapprox.(velocity.(cal,  X.R,  X.E, X.Rw, X.T, X.P, X.fluid), X.U))

    @test all(isapprox.(velocity.(cal,  X.R,  X.E, X.Rw, X.T, X.P, X.fluid),
                        velocity.(calt, X.R,  X.E, X.Rw, X.T, X.P, X.fluid)))

    @test all(isapprox.(velocity.(cal,  X.R,  X.E, X.Rw, X.T+10, X.P, X.fluid),
                        velocity.(calt, X.R,  X.E, X.Rw, X.T+10, X.P, X.fluid)))
    
    

    ΔT = X.Tw - cal.T

    k = heatcond(fluid, cal.T+theta*ΔT, cal.P)
    Pr = prandtl(fluid, cal.T+theta*ΔT, cal.P)
    nu = kinvisc(fluid, cal.T+theta*ΔT, cal.P)

    ϕ = X.Rw * k*Pr^n*ΔT
    
    
    A, B, n = king
    a0 = -A/B
    a1 =  1/B
    p = cal.fit.poly
    @test isapprox(a0, p[0]*nu^n)
    @test isapprox(a1, p[1]*nu^n/ϕ)
    
    

    # Test the corrections
    Em = mean(X.E) # Average voltage
    @test all(isapprox.(velocity.(cal,  X.R,  Em, X.Rw, X.T, X.P, X.fluid),
                        velocity.(calt, X.R,  Em, X.Rw, X.T, X.P, X.fluid)))

    @test all(isapprox.(velocity.(cal,  X.R,  Em, X.Rw*1.05, X.T+10,X.P*0.95,X.fluid),
                        velocity.(calt, X.R,  Em, X.Rw*1.05, X.T+10,X.P*0.95,X.fluid)))



    
end


# Test property changes

let
    n = 0.3
    theta = 0.5


    fluid  = ConstPropFluid(AIR, 293.15, 101325.0)
    # Let's change kin. visc without changing Pr, k
    fluid2 = ConstPropFluid(fluid.ρ, fluid.μ*2, fluid.k, fluid.cp/2) # Changed visc

    # Let's keep nu bu change k, Pr
    fluid3 = ConstPropFluid(fluid.ρ, fluid.μ, fluid.k*1.5, fluid.cp) # Changed visc and k
    
    X = example_calibr(fluid) #ConstPropFluid)
    king = X.king
    fit = E -> ( (E*E - king.A) / king.B )^(1/king.n)
    ifit = U -> king.A + king.B * U ^ king.n

    cal = HWCalibr(X.R, X.U, X.E, X.Rw, X.T, X.P, makekingfitfun(a=king.n, N=1),
                   fluid=fluid, n=n, theta=theta)

    
    ΔT = X.Tw - cal.T

    k = heatcond(fluid, cal.T+theta*ΔT, cal.P)
    Pr = prandtl(fluid, cal.T+theta*ΔT, cal.P)
    nu = kinvisc(fluid, cal.T+theta*ΔT, cal.P)

    ϕ = X.Rw * k*Pr^n*ΔT


    k2 = heatcond(fluid2, cal.T+theta*ΔT, cal.P)
    Pr2 = prandtl(fluid2, cal.T+theta*ΔT, cal.P)
    nu2 = kinvisc(fluid2, cal.T+theta*ΔT, cal.P)

    ϕ2 = X.Rw * k2*Pr2^n*ΔT

    k3 = heatcond(fluid3, cal.T+theta*ΔT, cal.P)
    Pr3 = prandtl(fluid3, cal.T+theta*ΔT, cal.P)
    nu3 = kinvisc(fluid3, cal.T+theta*ΔT, cal.P)

    ϕ3 = X.Rw * k3*Pr3^n*ΔT

    
    Em = mean(X.E) # Average voltage

    U0 = velocity(cal,  X.R,  Em, X.Rw, X.T, X.P, X.fluid)
    U2 = velocity(cal,  X.R,  Em, X.Rw, X.T, X.P, fluid2)

    @test isapprox(U2, U0 * (nu2/nu))

    f=1.05
    U0 = velocity(cal,  X.R,  Em, f*X.Rw, X.T+10, X.P, X.fluid)
    U2 = velocity(cal,  X.R,  Em, f*X.Rw, X.T+10, X.P, fluid2)

    @test isapprox(U2, U0*(nu2/nu)) 

    

    # Now we change k,Pr and keep nu constant

    U0 = velocity(cal,  X.R,  Em, X.Rw, X.T, X.P, X.fluid)
    U3 = velocity(cal,  X.R,  Em*sqrt(ϕ3/ϕ), X.Rw, X.T, X.P, fluid3)
    @test isapprox(U0, U3)

    # Now we change k,Pr and keep nu constant but change Rw and Tw
    U0 = velocity(cal,  X.R,  Em, X.Rw*1.05, X.T+10, X.P*0.95, X.fluid)
    U3 = velocity(cal,  X.R,  Em*sqrt(ϕ3/ϕ), X.Rw*1.05, X.T+10, X.P*0.95, fluid3)
    @test isapprox(U0, U3)
    
    
end


