using Statistics

# Testing wire

let D=5e-6, L=2e-3
    Ta = 273.15 + 20.0
    Pa = 101325.0
    AIRconst = ConstPropFluid(AIR, Ta, Pa)

    R = Resistor(R=3.5, a=0.4e-2, T=Ta)

    U = 1.0 : 10.0
    n = 0.42
    A = 2.0
    B = 1.0
    E = @. sqrt(A + B*U^n)
    
    fit = KingLaw(E, U, n)

    @test fit.A ≈ A
    @test fit.B ≈ B
    @test fit.n == n

    fit2 = KingLaw(E, U; n=0.45)
    @test fit2.A ≈ A
    @test fit2.B ≈ B
    @test fit2.n ≈ n
    
    Tw = 273.15 + 240.0
    Rw = R(Tw)

    calibr = TempCalibr(R, E, U, Ta, Pa, Rw, KingLaw; fluid=AIRconst)

    w = CTASensor(TempCalibr, R, Rw, E, U, Ta, Pa, KingLaw; fluid=AIRconst)
    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    

    U1 = w.(E)
    U2 = velocity(w, E)
    U3 = similar(E)
    velocity!(U3, w, E)
    
    @test all(U1 .≈ U)
    @test all(U1 .≈ U2)
    @test all(U2 .≈ U3)

    Em = mean(E)

    @test velocity(w, Em) ≈ velocity(calibr, Em)
    @test velocity(w, Em; T=Ta+10) ≈ velocity(calibr, Em; T=Ta+10)
    
    
end


let 
    Ta = 273.15 + 25.0
    Pa = 101325.0

    AIRconst = ConstPropFluid(AIR, Ta, Pa)

    R = Thermistor(R=5e3, B=3950.0, T=Ta)
    U = 1.0 : 10.0
    n = 0.35
    A = 36.0
    B = 10.0
    E = @. sqrt(A + B*U^n)
    
    fit = KingLaw(E, U, n)

    Rw = 100.0
    Tw = temperature(R, Rw)

    calibr = HWCalibr(R, E, U, Ta, Pa, Rw, (x,y)->PowerPoly(x,y; n=n, N=1);
                      fluid=AIRconst, theta=0.3, n=0.3)

    w =  CTASensor(calibr)

    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    
    U1 = w.(E)
    U2 = velocity(w, E)
    U3 = similar(E)
    velocity!(U3, w, E)
    @test all(U1 .≈ U)
    @test all(U1 .≈ U2)
    @test all(U2 .≈ U3)

    Em = mean(E)

    @test velocity(w, Em) ≈ velocity(calibr, Em)
    @test velocity(w, Em; T=Ta+10) ≈ velocity(calibr, Em; T=Ta+10)

end

