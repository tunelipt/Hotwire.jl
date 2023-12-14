
# Testing correction factor

let
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


let
    tc1 = TempCorrect(298.15, 100.0, 498.15)
    tc2 = TempCorrect(298.15+5, 100.0, 498.15)

    AIRconst = ConstPropFluid(AIR, 298.15, 101325.0)

    


end
