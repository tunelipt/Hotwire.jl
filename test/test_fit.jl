# Testing the curve fitting stuff

let
    U = 1.0:0.2:5.0
    N = 5
    a = 0.45
    E² = (2.0 .+ 2.0 * U.^a)


    f = KingPoly(E², U, 0.45, N)

    U1 = f.(E²)

    @test U1 ≈ U
    @test f.poly[0] ≈ -1.0
    @test f.poly[1] ≈ +0.5
    for i in 2:5
        @test f.poly[i] ≈ 0.0 atol=1e-10
    end
    

    mf = makekingfitfun(a=a, N=N)

    f2 = mf(E², U)
    
    U2 = f.(E²)

    @test U2 ≈ U
    @test f2.poly[0] ≈ -1.0
    @test f2.poly[1] ≈ +0.5
    for i in 2:5
        @test f2.poly[i] ≈ 0.0 atol=1e-10
    end
    
end
