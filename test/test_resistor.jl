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
    
