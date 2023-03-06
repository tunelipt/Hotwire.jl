# Testing Probe2d stuff

let
    trivialcal(E,T) = E
    wire1 = Wire(3.5, 6.3)
    wire2 = Wire(3.7, 6.7)

    probe = Probe2d((wire1, wire2), (trivialcal, trivialcal), (0.0, 0.0), "")

    Ux0 = 2.0; Uy0 = 0.0
    Uc1 = (Ux0 - Uy0)
    Uc2 = (Ux0 + Uy0)
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 atol=1e-5

    Ux0 = 2.0; Uy0 = 1.0
    Uc1 = (Ux0 - Uy0)
    Uc2 = (Ux0 + Uy0)
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0

    Ux0 = 2.0; Uy0 = -1.0
    Uc1 = (Ux0 - Uy0)
    Uc2 = (Ux0 + Uy0)
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0
    
    
end
