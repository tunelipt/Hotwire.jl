
# Testing correction factor


let
    # Testing temperature changes only
    tc1 = TempCorrect(298.15, 100.0, 498.15)
    tc2 = TempCorrect(298.15+5, 100.0, 498.15)
    tc3 = TempCorrect(298.15-5, 100.0, 498.15)
    tc4 = TempCorrect(298.15,  90.0, 498.15)
    
    @test anemcorrectfactor(tc1, tc1) ≈ 1.0
    @test anemcorrectfactor(tc2, tc2) ≈ 1.0
    @test anemcorrectfactor(tc3, tc3) ≈ 1.0
    
    
    @test anemcorrectfactor(tc1, tc2) ≈ sqrt( 200 / 195)
    @test anemcorrectfactor(tc1, tc3) ≈ sqrt( 200 / 205)
    @test anemcorrectfactor(tc1, tc4) ≈ sqrt( 100 / 90 )
    
    @test anemcorrectfactor(tc1, T=298.15+5, Rw=100.0, Tw=498.15) ≈ sqrt( 200 / 195)
    @test anemcorrectfactor(tc1, T=298.15-5, Rw=100.0, Tw=498.15) ≈ sqrt( 200 / 205)
    @test anemcorrectfactor(tc1, T=298.15, Rw=90.0, Tw=498.15) ≈ sqrt( 100 / 90 )
    
    
    @test tc1(tc2) ≈ sqrt( 200 / 195)
    @test tc1(tc3) ≈ sqrt( 200 / 205)
    @test tc1(tc4) ≈ sqrt( 100 / 90 )
    
    @test tc1(T=298.15+5, Rw=100.0, Tw=498.15) ≈ sqrt( 200 / 195)
    @test tc1(T=298.15-5, Rw=100.0, Tw=498.15) ≈ sqrt( 200 / 205)
    @test tc1(T=298.15, Rw=90.0, Tw=498.15) ≈ sqrt( 100 / 90 )
end


let # Testing TOnlyCorrect
    tc1 = TempCorrect(298.15, 100.0, 498.15)
    tc2 = TempCorrect(298.15+5, 100.0, 498.15)

    AIRconst = ConstPropFluid(AIR, 298.15, 101325.0)
    E = 1.0
    ν = kinvisc(AIRconst, 298.15, 101325.0)
    corr = TOnlyCorrect(ν)
    # Trivial case
    E1 = anemcorrect(E, tc1, corr, tc1, corr)
    @test E1 ≈ 1.0

    E1 = anemcorrect(E, tc1, corr, tc2, corr)
    @test E1 ≈ tc1(tc2)

    corr1 = TOnlyCorrect(corr, tc2, AIRconst, 101325.0)
    @test corr == corr1

end

    
let
    # Testing WireCorrect with constant properties
    T1 = 298.15
    R = Thermistor(R=5e3, B=3950.0, T=T1)
    T2 = T1 + 5.0
    P = 101325.0
    Rw = 100.0
    Tw = temperature(R, Rw)
    Rw1 = 90.0
    Tw1 = temperature(R, Rw1)
    
    tc1 = TempCorrect(T1, Rw, Tw)
    tc2 = TempCorrect(T2, Rw, Tw)
    tc3 = TempCorrect(T1, Rw1, Tw1)

    AIRconst = ConstPropFluid(AIR, T1, P)
    n=1/3
    corr1 = WireCorrect(tc1, AIRconst, P, n)
    corr2 = WireCorrect(corr1, tc2, AIRconst, P)
    corr3 = WireCorrect(corr1, tc3, AIRconst, P)
    
    E = 1.0
    @test anemcorrect(E, tc1, corr1, tc1, corr1) ≈ E
    @test anemcorrect(E, tc1, corr1, tc2, corr2) ≈ E*tc1(tc2)
    @test anemcorrect(E, tc1, corr1, tc3, corr3) ≈ E*tc1(tc3)
    
end


let
    # Testing WireCorrect with constant properties
    T1 = 298.15
    R = Thermistor(R=5e3, B=3950.0, T=T1)
    P = 101325.0
    Rw = 100.0
    Tw = temperature(R, Rw)
    Tm = (T1+Tw)/2
    tc1 = TempCorrect(T1, Rw, Tw)

    # Same temperatures but different fluids
    x1 = ConstPropFluid(AIR, T1, P)
    x2  = ConstPropFluid(HELIUM, T1, P)
    n = 1/3
    ϕ₁ = heatcond(x1, Tm, P) * prandtl(x1, Tm, P)^n
    ϕ₂ = heatcond(x2, Tm, P) * prandtl(x2, Tm, P)^n

    corr1 = WireCorrect(tc1, x1, P, n)
    corr2 = WireCorrect(tc1, x2, P, n)
    E = 1.0
    @test anemcorrect(E, tc1, corr1, tc1, corr2) ≈ E * sqrt(ϕ₁/ϕ₂)

end



let
    # Testing WireCorrect with variable properties
    T1 = 298.15
    T2 = T1 + 5.0
    P = 101325.0

    R = Thermistor(R=5e3, B=3950.0, T=T1)

    Rw = 100.0
    Tw = temperature(R, Rw)

    Rw1 = 90.0
    Tw1 = temperature(R, Rw1)

    Tm1 = (T1 + Tw) / 2
    Tm2 = (T2 + Tw) / 2
    Tm3 = (T1 + Tw1) / 2
    
    tc1 = TempCorrect(T1, Rw, Tw)
    tc2 = TempCorrect(T2, Rw, Tw)
    tc3 = TempCorrect(T1, Rw1, Tw1)

    n=1/3
    corr1 = WireCorrect(tc1, AIR, P, n)
    corr2 = WireCorrect(corr1, tc2, HELIUM, P)
    corr3 = WireCorrect(corr1, tc3, HELIUM, P)
    
    E = 1.0
    ϕ₁ = heatcond(AIR, Tm1, P) * prandtl(AIR, Tm1, P)^n
    ϕ₂ = heatcond(HELIUM, Tm2, P) * prandtl(HELIUM, Tm2, P)^n
    ϕ₃ = heatcond(HELIUM, Tm3, P) * prandtl(HELIUM, Tm3, P)^n

    @test ϕ₁ ≈ corr1.ϕ
    @test ϕ₂ ≈ corr2.ϕ
    @test ϕ₃ ≈ corr3.ϕ
    
    @test anemcorrect(E, tc1, corr1, tc1, corr1) ≈ E
    @test anemcorrect(E, tc1, corr1, tc2, corr2) ≈ E*sqrt(ϕ₁/ϕ₂)*tc1(tc2)
    @test anemcorrect(E, tc1, corr1, tc3, corr3) ≈ E*sqrt(ϕ₁/ϕ₃)*tc1(tc3)
    
end


let
    # Testing GlassbeadCorrect with constant properties
    T1 = 298.15
    T2 = T1 + 5
    P = 101325.0
    
    R = Thermistor(R=5e3, B=3950.0, T=T1)

    Rw = 100.0
    Tw = temperature(R, Rw)
    Tm1 = (T1+Tw)/2
    Tm2 = (T2+Tw)/2

    Rw1 = 90.0
    Tw1 = temperature(R, Rw1)
    Tm2 = (T1+Tw1)/2

    tc1 = TempCorrect(T1, Rw, Tw)
    tc2 = TempCorrect(T2, Rw, Tw)
    tc3 = TempCorrect(T1, Rw1, Tw1)

    # Same temperatures but different fluids
    x1 = ConstPropFluid(AIR, T1, P)
    n = 1/3

    corr1 = mf58correct(tc1, x1, P; n=n)
    corr2 = GlassbeadCorrect(corr1, tc2, x1, P)
    corr3 = GlassbeadCorrect(corr1, tc3, x1, P)
    
    E = 2.0
    @test anemcorrect(E, tc1, corr1, tc1, corr1) ≈ E
    @test anemcorrect(E, tc1, corr1, tc2, corr2) ≈ E*tc1(tc2)
    @test anemcorrect(E, tc1, corr1, tc3, corr3) ≈ E*tc1(tc3)
    

end


let
    # Testing Glass bead
    # If N=0 and β=0, it should be equal to wire corrrect
    T1 = 298.15
    T2 = T1 + 5.0
    P = 101325.0

    R = Thermistor(R=5e3, B=3950.0, T=T1)

    Rw = 100.0
    Tw = temperature(R, Rw)

    Rw1 = 90.0
    Tw1 = temperature(R, Rw1)

    Tm1 = (T1 + Tw) / 2
    Tm2 = (T2 + Tw) / 2
    Tm3 = (T1 + Tw1) / 2
    
    tc1 = TempCorrect(T1, Rw, Tw)
    tc2 = TempCorrect(T2, Rw, Tw)
    tc3 = TempCorrect(T1, Rw1, Tw1)

    n=1/3
    cw1 = WireCorrect(tc1, AIR, P, n)
    cw2 = WireCorrect(tc2, HELIUM, P, n)
    cw3 = WireCorrect(tc3, HELIUM, P, n)

    gb1 = mf58correct(tc1, AIR, P; n=n, N=0, beta=0.0)
    gb2 = GlassbeadCorrect(gb1, tc2, HELIUM, P)
    gb3 = GlassbeadCorrect(gb2, tc3, HELIUM, P)

    E = 2.0
    @test anemcorrect(E, tc1, gb1, tc1, gb1) ≈ E
    @test anemcorrect(E, tc1, gb1, tc2, gb2) ≈ anemcorrect(E, tc1, cw1, tc2, cw2)
    @test anemcorrect(E, tc1, gb1, tc3, gb3) ≈ anemcorrect(E, tc1, cw1, tc3, cw3)
    

end


let
    # Testing Glassbead: MF58, N=0, beta=150

    U = 0.5:0.1:10.0
    T1 = 298.15
    T2 = T1 + 5.0
    P = 101325.0

    R = Thermistor(R=5e3, B=3950.0, T=T1)

    Rw = 100.0
    Tw = temperature(R, Rw)

    Rw1 = 90.0
    Tw1 = temperature(R, Rw1)

    Tm1 = (T1 + Tw) / 2
    Tm2 = (T2 + Tw) / 2
    Tm3 = (T1 + Tw1) / 2
    
    tc1 = TempCorrect(T1, Rw, Tw)
    tc2 = TempCorrect(T2, Rw, Tw)
    tc3 = TempCorrect(T1, Rw1, Tw1)

    D = 2.0e-3
    L = 4.0e-3
    A = π*D*L
    
    N = 0
    n = 1/3

    ν₁ = kinvisc(AIR, Tm1, P)
    ν₂ = kinvisc(HELIUM, Tm2, P)
    ν₃ = kinvisc(HELIUM, Tm3, P)

    Re1 = U .* D ./ ν₁    
    Re2 = U .* D ./ ν₂    
    Re3 = U .* D ./ ν₃    

    Pr1 = prandtl(AIR, Tm1, P)
    Pr2 = prandtl(HELIUM, Tm2, P)
    Pr3 = prandtl(HELIUM, Tm3, P)

    Nu1 = 0.683 * Re1 .^ 0.466 .* Pr1 .^ n
    Nu2 = 0.683 * Re2 .^ 0.466 .* Pr2 .^ n
    Nu3 = 0.683 * Re3 .^ 0.466 .* Pr3 .^ n

    k1 = heatcond(AIR, Tm1, P)
    k2 = heatcond(HELIUM, Tm2, P)
    k3 = heatcond(HELIUM, Tm3, P)

    H1 = A * Nu1 .* k1 / D
    H2 = A * Nu2 .* k2 / D
    H3 = A * Nu3 .* k3 / D

    β = 150.0
    E1 = @. sqrt( H1 / (1 + β*H1)  * Rw * (Tw - T1) )
    E2 = @. sqrt( H2 / (1 + β*H2)  * Rw * (Tw - T2) )
    E3 = @. sqrt( H3 / (1 + β*H3)  * Rw1 * (Tw1 - T1) )

    lna, b  = Hotwire.linearfit(log.(U), ones(length(U)), log.(E))
    
    
    
    
    

end