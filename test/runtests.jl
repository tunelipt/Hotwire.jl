using Hotwire
using Polynomials
using Test

@testset "Hotwire.jl" begin

    # Resistance
    r = Resistor(1e3, 0.0, 298.15)

    @test resistance(r) == 1e3
    @test resistance(r, 298.15) == 1e3

    r = Resistor(1e3, 0.01, 298.15)
    @test resistance(r, 298.15) == 1e3
    @test resistance(r, 299.15) == 1.01e3

    @test temperature(r, 1010.0) == 299.15

    r = Resistor(1, 1//100, 0)
    @test resistance(r, 100) == 2//1
    @test temperature(r, 2) == 100//1
    
    # Thermistor
    r = Thermistor(1e3, 0.0, 298.15)
    @test resistance(r) == 1e3
    @test resistance(r, 303.15) == 1e3

    r = Thermistor(1e3, 3300.0, 298.15)
    @test temperature(r, resistance(r, 50.0+273.15)) ≈ 50.0+273.15
    @test temperature(r, resistance(r, 100.0+273.15)) ≈ 100.0+273.15
    @test 100.0 ≈ resistance(r, temperature(r, 100.0))

    r = Thermistor(5f3, 3.3f3, 298.15f0)
    @test temperature(r, resistance(r, 50f0+273.15f0)) ≈ 50f0+273.15f0
    @test temperature(r, resistance(r, 100f0+273.15f0)) ≈ 100f0+273.15f0
    @test 100f0 ≈ resistance(r, temperature(r, 100f0))
    


    # CTASensor

    r = Resistor(1.0, 0.01, 300.0)
    w = CTASensor(r, 2.0)
    @test r == resistor(w)
    @test optemperature(w) ≈ 400.0
    @test overtemp(w) ≈ 100.0
    @test overheat_ratio(w) ≈ 1.0

    w1 = Wire(1.0, 2.0, 300.0, 0.01)
    @test w == w1

    r = Resistor(1, 1//100, 0)
    w = CTASensor(r, 2//1)
    @test r == resistor(w)
    @test optemperature(w) == 100//1
    @test overtemp(w) ≈ 100//1
    @test overheat_ratio(w) ≈ 1//1
    
    
    r = Thermistor(1e3, 3300, 25.0+273.15)
    Rw = 100.0
    w = CTASensor(r, Rw)
    Tw = temperature(r, Rw)
    a = Rw/1e3 - 1.0
    @test r == resistor(w)
    @test optemperature(w) ≈ Tw
    @test overtemp(w) ≈ Tw - 298.15
    @test overheat_ratio(w) ≈ a

    r = Thermistor(1f3, 3300, 25f0+273.15f0)
    Rw = 100f0
    w = CTASensor(r, Rw)
    Tw = temperature(r, Rw)
    a = Rw/1f3 - 1f0
    @test r == resistor(w)
    @test optemperature(w) ≈ Tw
    @test overtemp(w) ≈ Tw - 298.15f0
    @test overheat_ratio(w) ≈ a
    
    # CCASensor - not much to do
    r = Resistor(1.0, 0.01, 300.0)
    w = CCASensor(r, 1.0)
    @test resistor(w) == r
    @test current(w) == 1.0

    # tempcorr.jl
    # tempcorr(CTASensor)
    Too = 20.0+273.15
    r = Resistor(3.5, 0.0036, Too)
    w = CTASensor(r, 5.5)
    @test 1.0 ≈ tempcorr(w, 1.0, Too) # No temperature correction
    @test tempcorr(w, 1.0, 25.0+273.15) > 1.0 # Increased temperature -> lower anemometer output
    @test tempcorr(w, 1.0, 15.0+273.15) < 1.0 # Decreased temperature -> higher anemometer output
    E1 = 1.0
    Too1 = 25.0+273.15
    E = tempcorr(w, E1, Too1)
    Tw = optemperature(w)
    @test E/E1 ≈ sqrt( (Tw - Too) / (Tw - Too1) )

    Too = 20f0+273.15f0
    r = Resistor(3.5f0, 0.0036f0, Too)
    w = CTASensor(r, 5.5f0)
    @test 1f0 ≈ tempcorr(w, 1.0f0, Too) # No temperature correction
    @test tempcorr(w, 1f0, 25f0+273.15f0) > 1f0 # Increased temperature -> lower anemometer output
    @test tempcorr(w, 1f0, 15f0+273.15f0) < 1f0 # Decreased temperature -> higher anemometer output
    E1 = 1f0
    Too1 = 25f0+273.15f0
    E = tempcorr(w, E1, Too1)
    Tw = optemperature(w)
    @test E/E1 ≈ sqrt( (Tw - Too) / (Tw - Too1) )
    
    Too = 20.0+273.15
    Too1 = 25.0+273.15
    r = Thermistor(1e3, 3300, Too)
    w = CTASensor(r, 10.0)
    @test 1.0 ≈ tempcorr(w, 1.0, 20.0+273.15) # No temperature correction
    @test tempcorr(w, 1.0, 25.0+273.15) > 1.0 # Increased temperature -> lower anemometer output
    @test tempcorr(w, 1.0, 15.0+273.15) < 1.0 # Decreased temperature -> higher anemometer output
    E1 = 1.0
    Too1 = 25.0+273.15
    E = tempcorr(w, E1, Too1)
    Tw = optemperature(w)
    @test E/E1 ≈ sqrt( (Tw - Too) / (Tw - Too1) )

    @test 1.0 ≈ tempcorr(w, 1.0, Too, Too)
    @test tempcorr(w, 1.0, 25.0+273.15, Too) > 1.0
    @test tempcorr(w, 1.0, 15.0+273.15, Too) < 1.0
    E = tempcorr(w, E1, Too1, Too)
    Tw = optemperature(w)
    @test E/E1 ≈ sqrt( (Tw - Too) / (Tw - Too1) )
    
    
    # tempcorr(CCASensor)
    To = 20.0+273.15
    r = Resistor(3.5, 0.0036, To)
    Two = 100.0+273.15

    Rwo = resistance(r, Two)
    I = 0.2
    Eo = Rwo * I
    hA = Eo*I / (Two - To)
    g = 1.0
    
    w = CCASensor(r, I, g)

    @test tempcorr(w, Eo, To) ≈ Eo
    @test tempcorr(w, Eo, To, To) ≈ Eo
    
    Te = 30.0+273.15
    E = Eo
    Eo = tempcorr(w, E, Te)
    Two = temperature(r, Eo/(g*I))
    Tw = temperature(r, E/(g*I))
    Qo = Eo/g * I
    Q = E/g * I
    hAo = Qo / (Two - To)
    hA = Q / (Tw - Te)
    @test Eo / E ≈ (Two - To) / (Tw - Te)
    @test hA ≈ hAo

    
    Te = 30.0+273.15
    E = Eo
    Eo = tempcorr(w, E, Te, temperature(w.R))
    Two = temperature(r, Eo/(g*I))
    Tw = temperature(r, E/(g*I))
    Qo = Eo/g * I
    Q = E/g * I
    hAo = Qo / (Two - To)
    hA = Q / (Tw - Te)
    @test Eo / E ≈ (Two - To) / (Tw - Te)
    @test hA ≈ hAo
    
    # Testing root finding equation
    @test Hotwire.funroot(x->(x^2-1.0), 0.5, 0.6) ≈ 1.0 atol=1e-4
    @test Hotwire.funroot(x->sin(x), π*1.1, π*1.02) ≈ π atol=1e-4

    # Testing tempcorr for CCASensor{Thermistor}
    To = 20.0+273.15
    r = Thermistor(1e3, 3300, To)
    g = 1.0
    I = 0.2
    w = CCASensor(r, I, g)
    @test tempcorr(w, 1.0, To) ≈ 1.0

    for Te in [15.0, 30.0] .+ 273.15
        for E in [0.5, 1.0, 2.0]
            Eo = tempcorr(w, E, Te)
            Two = temperature(r, Eo/(g*I))
            Tw = temperature(r, E/(g*I))
            @test Eo/E ≈ (Two-To) / (Tw - Te) atol=1e-6
            Qo = Eo/g * I
            Q = E/g * I
            hAo = Qo / (Two - To)
            hA = Q / (Tw - Te)
            @test hA ≈ hAo atol=1e-5
        end
    end

    for Te in [15.0, 30.0] .+ 273.15
        for E in [0.5, 1.0, 2.0]
            Eo = tempcorr(w, E, Te, To)
            Two = temperature(r, Eo/(g*I))
            Tw = temperature(r, E/(g*I))
            @test Eo/E ≈ (Two-To) / (Tw - Te) atol=1e-6
            Qo = Eo/g * I
            Q = E/g * I
            hAo = Qo / (Two - To)
            hA = Q / (Tw - Te)
            @test hA ≈ hAo atol=1e-5
        end
    end


    # tempcorr(CCASensor) Float32
    To = 20f0+273.15f0
    r = Resistor(3.5f0, 0.0036f0, To)
    Two = 100f0+273.15f0

    Rwo = resistance(r, Two)
    I = 0.2f0
    Eo = Rwo * I
    hA = Eo*I / (Two - To)
    g = 1f0
    
    w = CCASensor(r, I, g)

    @test tempcorr(w, Eo, To) ≈ Eo
    @test tempcorr(w, Eo, To, To) ≈ Eo
    
    Te = 30f0+273.15f0
    E = Eo
    Eo = tempcorr(w, E, Te)
    Two = temperature(r, Eo/(g*I))
    Tw = temperature(r, E/(g*I))
    Qo = Eo/g * I
    Q = E/g * I
    hAo = Qo / (Two - To)
    hA = Q / (Tw - Te)
    @test Eo / E ≈ (Two - To) / (Tw - Te)
    @test hA ≈ hAo

    
    Te = 30f0+273.15f0
    E = Eo
    Eo = tempcorr(w, E, Te, temperature(w.R))
    Two = temperature(r, Eo/(g*I))
    Tw = temperature(r, E/(g*I))
    Qo = Eo/g * I
    Q = E/g * I
    hAo = Qo / (Two - To)
    hA = Q / (Tw - Te)
    @test Eo / E ≈ (Two - To) / (Tw - Te)
    @test hA ≈ hAo
    
    # Testing root finding equation
    @test Hotwire.funroot(x->(x^2-1f0), 0.5f0, 0.6f0) ≈ 1.0 atol=1e-4
    @test Hotwire.funroot(x->sin(x), π*1.1f0, π*1.02f0) ≈ π atol=1e-4

    # Testing tempcorr for CCASensor{Thermistor}
    To = 20f0+273.15f0
    r = Thermistor(1f3, 3300, To)
    g = 1f0
    I = 0.2f0
    w = CCASensor(r, I, g)
    @test tempcorr(w, 1.0, To) ≈ 1f0

    for Te in [15f0, 30f0] .+ 273.15f0
        for E in [0.5f0, 1.0f0, 2.0f0]
            Eo = tempcorr(w, E, Te)
            Two = temperature(r, Eo/(g*I))
            Tw = temperature(r, E/(g*I))
            @test Eo/E ≈ (Two-To) / (Tw - Te) atol=1f-4
            Qo = Eo/g * I
            Q = E/g * I
            hAo = Qo / (Two - To)
            hA = Q / (Tw - Te)
            @test hA ≈ hAo atol=1f-4
        end
    end

    for Te in [15f0, 30f0] .+ 273.15f0
        for E in [0.5f0, 1.0f0, 2.0f0]
            Eo = tempcorr(w, E, Te, To)
            Two = temperature(r, Eo/(g*I))
            Tw = temperature(r, E/(g*I))
            @test Eo/E ≈ (Two-To) / (Tw - Te) atol=1f-4
            Qo = Eo/g * I
            Q = E/g * I
            hAo = Qo / (Two - To)
            hA = Q / (Tw - Te)
            @test hA ≈ hAo atol=1f-4
        end
    end

    

    # Testing calibration curve
    Rw = Thermistor(5e3, 3500, 25.0+273.15);
    cta = CTASensor(Rw, 100.0, 20/120);

    
    E = [1.333, 1.385, 1.423, 1.48, 1.521, 1.551, 1.575, 1.624]
    T0 = 20.0+273.15
    coefs = [2527.0, -3916.0, -1548.0, 6110.0, -3945.0, 815.0]
    fit = Polynomial(coefs)
    U = fit.(E)

    cal = calibr_curve(cta, U, E, T0, 5)
    
    @test cal.T0 == T0
    @test cal.fit.coeffs ≈ coefs rtol=1e-4
    @test all(cal.(cta, E) .≈ U)

    U = 1.0:0.5:10.0
    A = 2.0
    B = 0.7
    n = 0.42
    E = sqrt.(A .+ B .* U .^ n)
    fitk = KingLaw(E, U)
    
    cal = CalibrCurve(E, collect(U), fitk, T0)
    
    @test fitk.A ≈ A
    @test fitk.B ≈ B
    @test fitk.n ≈ n

    @test all(cal.(cta, E) .≈ U)
    
    
    include("test_probe2d.jl")
    include("test_probe3d.jl")
    
end
