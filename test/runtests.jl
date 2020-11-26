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
    w = CTASensor(r, 100.0)
    @test r == resistor(w)
    @test optemperature(w) == 100.0
    @test overtemp(w) == 100.0
    @test overheat_ratio(w) ≈ 1.0

    r = Thermistor(1e3, 3300, 25.0)
    w = CTASensor(r, 100.0)
    rw = resistance(r, 100.0)
    a = rw/1e3 - 1.0
    @test r == resistor(w)
    @test optemperature(w) == 100.0
    @test overtemp(w) ≈ 75.0
    @test overheat_ratio(w) ≈ a
    
    
end
