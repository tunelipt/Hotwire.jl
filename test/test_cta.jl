

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

    fit2 = KingLaw(E, U; n0=0.45)
    @test fit2.A ≈ A
    @test fit2.B ≈ B
    @test fit2.n ≈ n
    
    Tw = 273.15 + 240.0
    Rw = R(Tw)

    corr = TempCorrect(Ta, Pa, AIRconst, Rw, Tw)
    
    w =  CTASensor(R, Rw, Tw, 1.0, 0.0, corr, fit)

    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    
    ν_cal = kinvisc(corr)

    fc = CorrFactor(1.0, ν_cal)
    U1 = w.(E, fc)

    @test all(U1 .≈ U)

    # Let's change the temperatures and such

    tc1 = tempcorrect(Ta, Rw, Tw, Ta+5, Rw, Tw)
    E2 = E ./ tc1
    fc1 = correct.(w, E, T=Ta+5)

    # Since fluid properties are constant, ν e ν_cal are the same
    # and therefore the velocity is the same
    U2 = w.(E2, fc1)

    @test all(U2 .≈ U)

    Rw3 = Rw * 1.05
    Tw3 = temperature(R, Rw3)
    tc3 = tempcorrect(Ta, Rw, Tw, Ta+5, Rw3, Tw3)
    E3 = E ./ tc3
    fc3 = correct.(w, E, T=Ta+5, Rw=Rw3)
    U3 = w.(E3, fc3)

    @test all(U3 .≈ U)
    
    # Let's check if Probe1d interface works...
    probe = Probe1d(w, 0 #= setup can be anything =#)
    @test temperature(probe) == Tw
    @test resistance(probe) == Rw
    @test caltemp(probe) == Ta
    @test reftemp(probe) == Ta
    @test pressure(probe) == Pa

    U1b = probe.(E)
    @test all(U1b .≈ U)
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

    corr = mf58correct(Ta, Pa, AIRconst, Rw, Tw)    
    
    w =  CTASensor(R, Rw, Tw, 1.0, 0.0, corr, fit)

    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    
    ν_cal = kinvisc(corr)

    fc = CorrFactor(1.0, ν_cal)
    U1 = w.(E, fc)

    @test all(U1 .≈ U)

    # Let's change the temperatures and such

    tc1 = tempcorrect(Ta, Rw, Tw, Ta+5, Rw, Tw)
    E2 = E ./ tc1
    fc1 = correct.(w, E, T=Ta+5)

    # Since fluid properties are constant, ν e ν_cal are the same
    # and therefore the velocity is the same
    U2 = w.(E2, fc1)

    @test all(U2 .≈ U)

    Rw3 = Rw * 1.05
    Tw3 = temperature(R, Rw3)
    tc3 = tempcorrect(Ta, Rw, Tw, Ta+5, Rw3, Tw3)
    E3 = E ./ tc3
    fc3 = correct.(w, E, T=Ta+5, Rw=Rw3)
    U3 = w.(E3, fc3)

    @test all(U3 .≈ U)
    
    # Let's check if Probe1d interface works...
    probe = Probe1d(w, 0 #= setup can be anything =#)
    @test temperature(probe) == Tw
    @test resistance(probe) == Rw
    @test caltemp(probe) == Ta
    @test reftemp(probe) == Ta
    @test pressure(probe) == Pa

    U1b = probe.(E)
    @test all(U1b .≈ U)


end
