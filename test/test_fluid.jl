
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


let
    T = 298.15
    P = 101325.0
    fluid = "Helium"
    Helium = CPFluid(fluid)
    @test heatcond(Helium, T, P) == PropsSI("L", "T", T, "P", P, fluid)
    mu =  viscosity(Helium, T, P)
    rho = density(Helium, T, P)
    @test mu == PropsSI("V", "T", T, "P", P, fluid)
    @test kinvisc(Helium, T, P) == mu / rho
    @test specheat(Helium, T, P) == PropsSI("C", "T", T, "P", P, fluid)
    @test rho == PropsSI("D", "T", T, "P", P, fluid)
    @test prandtl(Helium, T, P) == PropsSI("PRANDTL", "T", T, "P", P, fluid)
    

end


let
    T = 298.15
    P = 101325.0
    R = 0.5

    fluid = HumidAir(T=T, P=P, R=R)
    
    w = HAPropsSI("W", "T", T, "P", P, "R", R)
    @test w == fluid.w
    
    
    @test heatcond(fluid, T, P) == PropsSI("L", "T", T, "P", P, "Air")
    mu =  viscosity(fluid, T, P)
    rho = density(fluid, T, P)
    @test mu == PropsSI("V", "T", T, "P", P, "Air")
    @test PropsSI("V", "T", T, "P", P, "Air") / rho == kinvisc(fluid, T, P)
    @test specheat(fluid, T, P) == HAPropsSI("cp_ha", "T", T, "P", P, "R", R)
    @test rho == 1/HAPropsSI("Vha", "T", T, "P", P, "R", R)
    @test prandtl(fluid, T, P) == PropsSI("PRANDTL", "T", T, "P", P, "Air")
    

end
