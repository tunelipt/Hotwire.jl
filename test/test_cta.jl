using Statistics

# Testing wire

let D=5e-6, L=2e-3
    Ta = 273.15 + 20.0
    Pa = 101325.0
    signal = linsignal(1,0)
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
    
    w =  CTASensor(R, Rw, Tw, signal, corr, fit)

    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    
    ν_cal = kinvisc(corr)

    fc = CorrFactor(1.0, ν_cal, E)
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
    signal = linsignal(1,0)

    w =  CTASensor(R, Rw, Tw, signal, corr, fit)

    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    
    ν_cal = kinvisc(corr)

    fc = CorrFactor(1.0, ν_cal, E)
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


# Testing Probe1d
let
    Ta = 273.15 + 20.0
    Pa = 101325.0
    signal = linsignal(1,0)
    AIRconst = ConstPropFluid(AIR, Ta, Pa)
    R = Resistor(R=3.5, a=0.4e-2, T=Ta)

    U = 1.0 : 10.0
    n = 0.42
    A = 2.0
    B = 1.0
    E = @. sqrt(A + B * U ^ n)
    ufun(E) = ( (E*E - A) / B )^(1/n)
    
    # Overheat ratio
    a = 0.8
    Rw = (1+a) * refresist(R)
    Tw = temperature(R, Rw)
    corr = TempCorrect(Ta, Pa, AIRconst, Rw, Tw)
    cta =  CTASensor(R, Rw, Tw, signal, corr, ufun)
    w = Probe1d(cta, (1,2,3))
    
    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test pressure(w) == Pa
    @test resistor(w) == R
    @test gain(w) == 1
    @test offset(w) == 0
    @test overheatratio(w) ≈ a
    @test sensorvolt(w, 1) == 1
    @test outsignal(w, 1) == 1
    
    ν_cal = kinvisc(corr)

    fc = CorrFactor(1.0, ν_cal, E)
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

    @test all(w.(E2, T=Ta+5) .≈ U)
    @test all(w.(E2, T=Ta+5, P=Pa*1.05) .≈ U)  # Constant properties, shouldn't change
    @test all(w.(E2, T=Ta+5, P=Pa*1.05, Rw=Rw, fluid=AIRconst)
              .≈ U)  # Constant properties, shouldn't change
    @test all(velocity.(w, E2, T=Ta+5, P=Pa*1.05, Rw=Rw, fluid=AIRconst)
              .≈ U)  # Constant properties, shouldn't change
    
    # Let's change the overheat:
    a2 = 0.85
    Rw2 = refresist(R) * (1 + a2)
    Tw2 = temperature(R, Rw2)
    tc2 = tempcorrect(Ta, Rw, Tw, Ta+5, Rw2, Tw2)
    E3 = E ./ tc2
    fc2 = correct.(w, E, T=Ta+5, Rw=Rw2)
    @test all(w.(E3, fc2) .≈ U)

    @test all(velocity.(w, E3, T=Ta+5, P=Pa*1.05, Rw=Rw2, fluid=AIRconst)
              .≈ U)  # Constant properties, shouldn't change
    
    
end


# Testing dantec1d
let
    Ta = 20.0
    Tk = 273.15 + Ta
    Pa = 101325.0
    signal = linsignal(1,0)
    AIRconst = ConstPropFluid(AIR, Tk, Pa)
    R = Resistor(R=3.5, a=0.4e-2, T=Tk)

    U = 1.0 : 10.0
    n = 0.42
    A = 2.0
    B = 1.0
    E = @. sqrt(A + B * U ^ n)
    ufun(E) = ( (E*E - A) / B )^(1/n)
    Nu = length(U)
    # Calibration table:
    caltable = [U E fill(Ta, Nu) fill(Pa/1000, Nu) E U]

    # Hardware config
    hconfig = [3.0, 0.8, 151.772, 4.528, 0.2, 3.828, Ta, Ta, 0.0, 1.0, 0.0, 0.0, 1.0,
               3.0, 8.0, 0.0]

    fitfun(E,U) = KingLaw(E,U,n)  # Should be exact...
    
    w = dantec1d(hconfig, caltable, fitfun; fluid=AIRconst)

    @test overheatratio(w) .≈ 0.8 atol=0.002  # The table is rounded...

    # Let's test the differing methods calling stuff
    cta = sensor(w)

    @test all(velf.(cta, E) .≈ U)
    @test all(velf.(w, E) .≈ U)
    @test all([velocity(cta, e, T=Tk, P=Pa, fluid=AIRconst, Rw=resistance(cta))
               for e in E] .≈ U)
    @test all([velocity(w, e, T=Tk, P=Pa, fluid=AIRconst, Rw=resistance(cta))
               for e in E] .≈ U)

    fc = correct(w, E[end])
    @test fc.f ≈ 1
    @test fc.nu == kinvisc(cta.corr)
    @test fc.E ≈ E[end]

    E1 = fill(E[end], 20)
    U1 = U[end]

    @test all(velocity.(cta, E1, fc) .≈ U1)
    @test all(velocity.(w, E1, fc) .≈ U1)

    Ux = zeros(length(E1))
    velocity!(Ux, w, E1)
    @test all(Ux .≈ U1)

    Rw = resistance(w)
    Tw = temperature(w)
    T0  = reftemp(w)
    R0 = refresist(w)
    R = resistor(w)
    
    @test Rw == resistance(cta)
    @test Tw == temperature(cta)
    @test T0  == reftemp(cta)
    @test R0  == refresist(cta)
    @test R  == resistor(cta)

    a2 = 0.85
    Rw2 = R0*(1+a2)
    Tw2 = temperature(R, Rw2)
    T2 = T0 + 5
    tc1 = tempcorrect(T0, Rw, Tw, T2, Rw2, Tw2)
    E2 = E1 ./ tc1
    
    velocity!(Ux, w, E2; T=T2, Rw=Rw2)
    @test all(Ux .≈ U1)
    @test all(velocity(w, E2, T=T2, Rw=Rw2) .≈ U1)
    @test all(w.(E2, T=T2, Rw=Rw2) .≈ U1)

    fc = correct(cta, mean(E2), T=T2, Rw=Rw2)
    @test all(velocity.(w, E2, fc) .≈ U1)
    @test all(w.(E2, fc) .≈ U1)
    
    
    
    
end

