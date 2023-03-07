# Testing Probe2d stuff

let
    E = collect(1.0:20.0)
    trivialcal = CalibrCurve(E, E, x->x, 293.15)
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

    k1 = 0.04; k2 = 0.04
    probe = Probe2d((wire1, wire2), (trivialcal, trivialcal), (k1, k2), "")

    Ux0 = 2.0; Uy0 = 0.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 atol=1e-5

    Ux0 = 2.0; Uy0 = 1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 

    Ux0 = 2.0; Uy0 = -1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 
    
    
    k1 = 0.03; k2 = 0.05
    probe = Probe2d((wire1, wire2), (trivialcal, trivialcal), (k1, k2), "")

    Ux0 = 2.0; Uy0 = 0.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 atol=1e-5

    Ux0 = 2.0; Uy0 = 1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 

    Ux0 = 2.0; Uy0 = -1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, 293.15)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 

    k1 = 0.03; k2 = 0.05
    probe = Probe2d((wire1, wire2), (trivialcal, trivialcal), (k1, k2), "")

    ang = -40.0:5.0:40.0
    Uc = 10.0
    Ux0 = Uc .* cosd.(ang)
    Uy0 = Uc .* sind.(ang)
    U1 = @. sqrt(2)/2 * (Ux0 + Uy0)
    U2 = @. sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = @. sqrt(k1*U1^2 +    U2^2)
    Ue2 = @. sqrt(   U1^2 + k2*U2^2)
    Uc1 = @. sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = @. sqrt( 2*Ue2^2 / (1 + k2))

    k1b, k2b = dircalibr(probe, ang, Uc, Uc1, Uc2)

    @test k1b ≈ k1
    @test k2b ≈ k2
    
end
