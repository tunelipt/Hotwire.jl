
let
    AIRconst = ConstPropFluid(AIR, 298.15, 101325.0)
    T = 298.15
    P = 101325.0
    @test heatcond(AIRconst, T, P) == heatcond(AIR, T, P)
    @test viscosity(AIRconst, T, P) == viscosity(AIR, T, P)
    @test kinvisc(AIRconst, T, P) == kinvisc(AIR, T, P)
    @test specheat(AIRconst, T, P) == specheat(AIR, T, P)
    @test density(AIRconst, T, P) == density(AIR, T, P)
    @test prandtl(AIRconst, T, P) == prandtl(AIR, T, P)
end
