# Testing probe3d stuff

let
    Ta = 273.15 + 20.0
    Pa = 101325.0

    AIRconst = ConstPropFluid(AIR, Ta, Pa)

    k² = [0.0, 0.0, 0.0]
    h² = [1.0, 1.0, 1.0]

    R = Resistor(R=3.5, a=0.4e-2, T=Ta)
    E = collect(1.0:20.0)
    Tw = 273.15 + 240.0
    Rw = R(Tw)

    corr = TempCorrect(Ta, Pa, AIRconst, Rw, Tw)

    wire1 = CTASensor(R, Rw, Tw, 1.0, 0.0, corr, E->E)
    wire2 = CTASensor(R, Rw, Tw, 1.0, 0.0, corr, E->E)
    wire3 = CTASensor(R, Rw, Tw, 1.0, 0.0, corr, E->E)
    wires = (wire1, wire2, wire3)

    probe = Probe3d(wires, k², h²)
    ν = kinvisc(corr)
    fc = CorrFactor((1.0, 1.0, 1.0), (ν, ν, ν), (1.0, 1.0, 1.0))

    Ux, Uy, Uz = probe(1.0, 1.0, 1.0, fc)
    @test Ux ≈ 1.0
    @test Uy ≈ 0.0 atol=1e-5
    @test Uz ≈ 0.0 atol=1e-5
    
    k² = [0.03, 0.04, 0.05]
    h² = [1.05, 1.08, 1.02]
    probe = Probe3d(wires, k², h²)

    Ux, Uy, Uz = probe(1.0, 1.0, 1.0, fc)
    @test Ux ≈ 1.0
    @test Uy ≈ 0.0 atol=1e-5
    @test Uz ≈ 0.0 atol=1e-5

    # Let's test directional calibration
    k² = [0.03, 0.04, 0.05]
    h² = [1.05, 1.08, 1.02]
    probe = Probe3d(wires, k², h²)

    
    θ = 30.0
    ϕ = 0.0:10.0:350.0
    Uc = 10.0
    Ux = Uc .* cosd(θ) .+ 0*ϕ
    uperp = Uc * sind(θ)
    Uz,Uy = uperp * sind.(ϕ), -uperp * cosd.(ϕ)
    cosw = probe.cosϕ
    U₁² = (Ux * cosw[1,1] + Uy * cosw[1,2] + Uz * cosw[1,3]).^2
    U₂² = (Ux * cosw[2,1] + Uy * cosw[2,2] + Uz * cosw[2,3]).^2
    U₃² = (Ux * cosw[3,1] + Uy * cosw[3,2] + Uz * cosw[3,3]).^2
    
    K = ones(3,3)
    for i in 1:3
        K[i,i] = probe.k²[i]
        K[i,probe.ih[i]] = probe.h²[i]
    end


    
    Ue1² = K[1,1]*U₁² +  K[1,2]*U₂² + K[1,3]*U₃² 
    Ue2² = K[2,1]*U₁² +  K[2,2]*U₂² + K[2,3]*U₃² 
    Ue3² = K[3,1]*U₁² +  K[3,2]*U₂² + K[3,3]*U₃²
    c = probe.c2e
    Uc1 = sqrt.( Ue1² ./ c[1] )
    Uc2 = sqrt.( Ue2² ./ c[2] )
    Uc3 = sqrt.( Ue3² ./ c[3] )
    
                
    vel1 = probe.(Uc1, Uc2, Uc3, fc)
    Ux1 = [v[1] for v in vel1]
    Uy1 = [v[2] for v in vel1]
    Uz1 = [v[3] for v in vel1]

    @test Ux1 ≈ Ux
    @test Uy1 ≈ Uy
    @test Uz1 ≈ Uz

    k₁², h₁² = dircalibr(probe, Uc, Uc1, Uc2, Uc3, ϕ, θ)

    @test k₁² ≈ k²
    @test h₁² ≈ h²
    
end
