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
    
    
    Tw = 273.15 + 240.0
    Rw = R(Tw)
    ΔT = Tw - Ta
    
    calibr = TempCalibr(R, U, E, Rw, Ta, Pa, makekingfitfun(a=n, N=1); fluid=AIRconst)

      
    sensor = CTASensor(R, Rw, calibr)

    w = Probe1d(sensor, "Cofiguration of the anememoter", ("dev1", 1))
    
    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test calpress(w) == Pa
    @test resistor(w) == R
    

    U1 = w.(E)
    U2 = velocity(w, E)
    U3 = similar(E)
    velocity!(U3, w, E)
    
    @test all(U1 .≈ U)
    @test all(U1 .≈ U2)
    @test all(U2 .≈ U3)

    Em = mean(E)

    @test velocity(w, Em) ≈ velocity(calibr, R, Em, Rw, Ta, Pa, AIRconst)
    @test velocity(w, Em; T=Ta+10) ≈ velocity(calibr, R, Em, Rw, Ta+10, Pa, AIRconst)
    
    
end


function xlet()
    Ta = 273.15 + 25.0
    Pa = 101325.0
    
    AIRconst = AIR #ConstPropFluid(AIR, Ta, Pa)
    fl = AIRconst
    
    R = Thermistor(R=5e3, B=3950.0, T=Ta)
    U = 1.0 : 10.0
    n = 0.4
    A = 36.0
    B = 10.0
    E = @. sqrt(A + B*U^n)
    

    Rw = 100.0
    Tw = temperature(R, Rw)
    ΔT = Tw - Ta

    theta=0.3
    n = 0.4
    Tf = Ta + theta*ΔT

    nu = kinvisc(AIRconst, Tf, Pa)
    k  = heatcond(AIRconst, Tf, Pa)
    Pr = prandtl(AIRconst, Tf, Pa)
    phi = k*Pr^n*Rw*ΔT
    
    calibr = HWCalibr(R, U, E, Rw, Ta, Pa, makekingfitfun(a=n, N=1);
                      fluid=AIRconst, n=n, theta=theta)

    p = calibr.fit.poly
    @test p[0] ≈ -A/(B*nu^n)
    @test p[1]/phi ≈ 1/(B*nu^n)
    
    
    w =  CTASensor(R, Rw, calibr)

    
    @test temperature(w) == Tw
    @test resistance(w) == Rw
    @test caltemp(w) == Ta
    @test calpress(w) == Pa
    @test resistor(w) == R
    
    U1 = w.(E)
    U2 = velocity(w, E)
    U3 = similar(E)
    velocity!(U3, w, E)
    @test all(U1 .≈ U)
    @test all(U1 .≈ U2)
    @test all(U2 .≈ U3)

    Em = mean(E)

    @test velocity(w, Em) ≈ velocity(calibr, R, Em, Rw, Ta, Pa, AIRconst)
    @test velocity(w, Em; T=Ta+10) ≈ velocity(calibr, R, Em, Rw, Ta+10, Pa, AIRconst)

end

