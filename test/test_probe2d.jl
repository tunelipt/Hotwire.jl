# Testing Probe2d stuff

let
    Ta = 273.15 + 20.0
    Pa = 101325.0

    AIRconst = ConstPropFluid(AIR, Ta, Pa)

    R = Resistor(R=3.5, a=0.4e-2, T=Ta)
    E = collect(1.0:20.0)
    Tw = 273.15 + 240.0
    Rw = R(Tw)

    corr = TempCorrect(Ta, Pa, AIRconst, Rw, Tw)
    ν = 
    wire1 = CTASensor(R, Rw, Tw, 1.0, corr, E->E)
    wire2 = CTASensor(R, Rw, Tw, 1.0, corr, E->E)
    

    probe = Probe2d((wire1, wire2), (0.0, 0.0), "")

    Ux0 = 2.0; Uy0 = 0.0
    Uc1 = (Ux0 - Uy0)
    Uc2 = (Ux0 + Uy0)
    fc = CorrFactor((1.0, 1.0), (kinvisc(corr), kinvisc(corr)))

    # E1 = Uc1, E2 = Uc2 - simplifiying...
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 atol=1e-5

    Ux0 = 2.0; Uy0 = 1.0
    Uc1 = (Ux0 - Uy0)
    Uc2 = (Ux0 + Uy0)
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0

    Ux0 = 2.0; Uy0 = -1.0
    Uc1 = (Ux0 - Uy0)
    Uc2 = (Ux0 + Uy0)
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0

    k1 = 0.04; k2 = 0.04
    probe = Probe2d((wire1, wire2), (k1, k2), "")

    Ux0 = 2.0; Uy0 = 0.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 atol=1e-5

    Ux0 = 2.0; Uy0 = 1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 

    Ux0 = 2.0; Uy0 = -1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 
    
    
    k1 = 0.03; k2 = 0.05
    probe = Probe2d((wire1, wire2), (k1, k2), "")

    Ux0 = 2.0; Uy0 = 0.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 atol=1e-5

    Ux0 = 2.0; Uy0 = 1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 

    Ux0 = 2.0; Uy0 = -1.0
    U1 = sqrt(2)/2 * (Ux0 + Uy0)
    U2 = sqrt(2)/2 * (Ux0 - Uy0)
    Ue1 = sqrt(k1*U1^2 +    U2^2)
    Ue2 = sqrt(   U1^2 + k2*U2^2)
    Uc1 = sqrt( 2*Ue1^2 / (1 + k1))
    Uc2 = sqrt( 2*Ue2^2 / (1 + k2))
    Ux, Uy = velocity(probe, Uc1, Uc2, fc)
    @test Ux ≈ Ux0
    @test Uy ≈ Uy0 
    
    k1 = 0.03; k2 = 0.05

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

    k1b, k2b = dircalibr(ang, Uc, Uc1, Uc2)

    @test k1b ≈ k1
    @test k2b ≈ k2
    
end
