
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
    
    return (E=E, U=U, T=Tc, P=Pc, Rw=Rw, fluid=fluid, R=R, king=king)
end


let
    X = example_calibr(AIR) #ConstPropFluid)
    king = X.king
    # Let's test King's law curve fit
    fit = KingLaw(X.E, X.U)

    @test isapprox(fit.A, king.A, atol=1e-7)
    @test isapprox(fit.B, king.B, atol=1e-7)
    @test isapprox(fit.n, king.n, atol=1e-7)

    Tm, Pm, caltable = Hotwire.makecaltable(X.E, X.U, X.T, X.P)
    Tw = temperature(X.R, X.Rw)
    cal0 = TempCalibr(X.R, X.Rw, Tw, X.T, X.P, X.fluid,
                      caltable, fit)
    @test isapprox(X.T, cal0.T)
    @test isapprox(X.P, cal0.P)
    @test isapprox(X.Rw, cal0.Rw)
    @test isapprox(Tw, cal0.Tw)
    @test isapprox(fit.A, cal0.fit.A)
    @test isapprox(fit.B, cal0.fit.B)
    @test isapprox(fit.n, cal0.fit.n)
    
    cal1 = TempCalibr(X.R, X.E, X.U, X.T, X.P, X.Rw, KingLaw)
    @test isapprox(X.T, cal1.T)
    @test isapprox(X.P, cal1.P)
    @test isapprox(X.Rw, cal1.Rw)
    @test isapprox(Tw, cal1.Tw)
    @test isapprox(fit.A, cal1.fit.A)
    @test isapprox(fit.B, cal1.fit.B)
    @test isapprox(fit.n, cal1.fit.n)

    N = length(X.E)
    cal = TempCalibr(X.R, X.E, X.U, ones(N)*X.T, ones(N)*X.P, ones(N)*X.Rw, KingLaw)
    @test isapprox(cal.T, cal1.T)
    @test isapprox(cal.P, cal1.P)
    @test isapprox(cal.Rw, cal1.Rw)

    # Test variations in Operating resistance
    # We will create a new calibration where Rw is changed by a fator of f.
    f = 1.1
    calx = TempCalibr(X.R, X.E, X.U, X.T, X.P, X.Rw, KingLaw; Rw=cal.Rw*f)
    α = X.R.α
    @test isapprox(calx.Tw-calx.T, (f-1)/α + f * (cal.Tw - cal.T))

    @test all(isapprox.(velocity.(cal, X.E), velocity.(calx, X.E; Rw=X.Rw)))
    @test all(isapprox.(velocity.(cal, X.E; Rw=X.Rw*f), velocity.(calx, X.E)))

    # Changes in pressure: we should not have any!
    @test all(isapprox.(velocity.(cal, X.E), velocity.(cal, X.E, P=X.P*0.9)))

    # Check changes in temperature
    Ta = cal.T + 10.0
    E1 = X.E .* sqrt.( (cal.Tw - Ta) / (cal.Tw - cal.T) )
    U1 = velocity.(cal, E1; T=Ta)
    U2 = velocity.(cal, X.E)
    @test all(isapprox.(U1, U2))

    N = length(X.E)
    Ta = 293.15 + 10.0
    U1 = velocity.(cal, X.E; T=Ta)
    U2 = velocity(cal, X.E; T=Ta)
    U3 = zeros(N)
    velocity!(U3, cal, X.E; T=Ta)

    @test all(U1 .== U2)
    @test all(U1 .== U3)
    @test all(U2 .== U3)
    
    
end


# Now we will test HWCalibr - using Constant properties!
let
    n = 0.3
    theta = 0.5
    
    fluid  = ConstPropFluid(AIR, 293.15, 101325.0)
    #luid1 = Co
    #luid2 = ConstPropFluid(CO2, 293.15, 101325.0)

    X = example_calibr(fluid) #ConstPropFluid)
    king = X.king
    # Let's test King's law curve fit
    fit = KingLaw(X.E, X.U)

    # Correct only temperature. Since properties are constant,
    # it should be identical to TempCalibr
    calt = TempCalibr(X.R, X.E, X.U, X.T, X.P, X.Rw, KingLaw; fluid=fluid)
        
    cal = HWCalibr(X.R, X.E, X.U, X.T, X.P, X.Rw, KingLaw;
                   fluid=fluid, n=n, theta=theta)
    @test isapprox(X.T, cal.T)
    @test isapprox(X.P, cal.P)
    @test isapprox(X.Rw, cal.Rw)
    @test isapprox(Tw, cal.Tw)

    N = length(X.E)
    cal = TempCalibr(X.R, X.E, X.U, ones(N)*X.T, ones(N)*X.P, ones(N)*X.Rw, KingLaw)
    @test isapprox(cal.T, cal1.T)
    @test isapprox(cal.P, cal1.P)
    @test isapprox(cal.Rw, cal1.Rw)

    # Test variations in Operating resistance
    # We will create a new calibration where Rw is changed by a fator of f.
    f = 1.1
    calx = TempCalibr(X.R, X.E, X.U, X.T, X.P, X.Rw, KingLaw; Rw=cal.Rw*f)
    α = X.R.α
    @test isapprox(calx.Tw-calx.T, (f-1)/α + f * (cal.Tw - cal.T))

    @test all(isapprox.(velocity.(cal, X.E), velocity.(calx, X.E; Rw=X.Rw)))
    @test all(isapprox.(velocity.(cal, X.E; Rw=X.Rw*f), velocity.(calx, X.E)))

    # Changes in pressure: we should not have any!
    @test all(isapprox.(velocity.(cal, X.E), velocity.(cal, X.E, P=X.P*0.9)))

    # Check changes in temperature
    Ta = cal.T + 10.0
    E1 = X.E .* sqrt.( (cal.Tw - Ta) / (cal.Tw - cal.T) )
    U1 = velocity.(cal, E1; T=Ta)
    U2 = velocity.(cal, X.E)
    @test all(isapprox.(U1, U2))
end
