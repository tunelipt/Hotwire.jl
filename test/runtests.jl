using Hotwire
using Test

@testset "Hotwire.jl" begin

    # Resistance
    r = Resistor(1e3, 0.0, 20.0)

    @test resistance(r) == 1e3
    @test resistance(r, 20.0) == 1e3

    r = Resistor(1e3, 1.0, 20.0)
    @test resistance(r, 20.0) == 1e3
    @test resistance(r, 21.0) == 1.01e3

    @test temperature(r, 1010.0) == 21.0

    # Thermistor
    r = Thermistor(1e3, 0.0, 25.0)
    @test resistance(r) == 1e3
    @test resistance(r, 30.0) == 1e3

    r = Thermistor(1e3, 3300.0, 25.0)
    @test temperature(r, resistance(r, 50.0)) ≈ 50.0
    @test temperature(r, resistance(r, 100.0)) ≈ 100.0

    @test 100.0 ≈ resistance(r, temperature(r, 100.0))


    # CTASensor

    r = Resistor(1.0, 1.0, 0.0)
    w = CTASensor(r, 2.0)
    @test r == resistor(w)
    @test optemperature(w) ≈ 100.0
    @test overtemp(w) ≈ 100.0
    @test overheat_ratio(w) ≈ 1.0

    r = Thermistor(1e3, 3300, 25.0)
    Rw = 100.0
    w = CTASensor(r, Rw)
    Tw = temperature(r, Rw)
    a = Rw/1e3 - 1.0
    @test r == resistor(w)
    @test optemperature(w) ≈ Tw
    @test overtemp(w) ≈ Tw - 25.0
    @test overheat_ratio(w) ≈ a

    # CCASensor - not much to do
    r = Resistor(1.0, 1.0, 0.0)
    w = CCASensor(r, 1.0)
    @test resistor(w) == r
    @test current(w) == 1.0

    # tempcorr.jl
    # tempcorr(CTASensor)
    Too = 20.0
    r = Resistor(3.5, 0.36, Too)
    w = CTASensor(r, 5.5)
    @test 1.0 ≈ tempcorr(w, 1.0, Too) # No temperature correction
    @test tempcorr(w, 1.0, 25.0) > 1.0 # Increased temperature -> lower anemometer output
    @test tempcorr(w, 1.0, 15.0) < 1.0 # Decreased temperature -> higher anemometer output
    E1 = 1.0
    Too1 = 25.0
    E = tempcorr(w, E1, Too1)
    Tw = optemperature(w)
    @test E/E1 ≈ sqrt( (Tw - Too) / (Tw - Too1) )
    
    
    r = Thermistor(1e3, 3300, Too)
    w = CTASensor(r, 10.0)
    @test 1.0 ≈ tempcorr(w, 1.0, 20.0) # No temperature correction
    @test tempcorr(w, 1.0, 25.0) > 1.0 # Increased temperature -> lower anemometer output
    @test tempcorr(w, 1.0, 15.0) < 1.0 # Decreased temperature -> higher anemometer output
    E1 = 1.0
    Too1 = 25.0
    E = tempcorr(w, E1, Too1)
    Tw = optemperature(w)
    @test E/E1 ≈ sqrt( (Tw - Too) / (Tw - Too1) )

    # tempcorr(CCASensor)
    To = 20.0
    r = Resistor(3.5, 0.36, To)
    Two = 100.0

    Rwo = resistance(r, Two)
    I = 0.2
    Eo = Rwo * I
    hA = Eo*I / (Two - To)
    g = 1.0
    
    w = CCASensor(r, I, g)

    @test tempcorr(w, Eo, To) ≈ Eo

    Te = 30.0
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


    # Testing root finding equation
    @test Hotwire.funroot(x->(x^2-1.0), 0.5, 0.6) ≈ 1.0 atol=1e-4
    @test Hotwire.funroot(x->sin(x), π*1.1, π*1.02) ≈ π atol=1e-4

    # Testing tempcorr for CCASensor{Thermistor}
    To = 20.0
    r = Thermistor(1e3, 3300, To)
    g = 1.0
    I = 0.2
    w = CCASensor(r, I, g)
    @test tempcorr(w, 1.0, To) ≈ 1.0

    for Te in [15.0, 30.0]
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
    
    
end
