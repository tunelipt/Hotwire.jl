using Hotwire
using Polynomials
using Test

function nucyl(Re,Pr)
    if Re < 4
        CRem = 0.989 * Re^0.33
    elseif Re < 40
        CRem = 0.911 * Re^0.385
    elseif Re < 4000
        CRem = 0.683 * Re^0.466
    elseif Re < 40_000
        CRem = 0.193 * Re^0.618
    else
        CRem = 0.027 * Re^0.805
    end
    return CRem*cbrt(Pr)
end
    
@testset "Hotwire.jl" begin

    include("test_resistor.jl")
    include("test_fluid.jl")
    include("test_calibr.jl")
    include("test_cta.jl")
    include("test_probe2d.jl")
    #include("test_probe3d.jl")
    
    
end
