using Polynomials


    
struct KingPoly{T}
    poly::Polynomial{T}
    a::T
end

Base.broadcastable(fit::KingPoly) = Ref(fit)

"""
`KingPoly(x, y, α=0.45; N=5)`

Fits a King Polynomial to the calibration curve of a thermal anemometer.

King's law has strong connection to the physics of the problem bu it
is often not accurate enough. On the other hand, a simple polynomial
between `E` and `U` does work well but has all the problems associated
with polynomials. A possibility is to get the best of both worlds:
expand King's law!

The proposed approximation is:

```math
Uᵅ = a₀ + a₁⋅E² + a₂⋅(E²)² + a₃⋅(E²)³ + ⋯ + aₚ⋅(E²)ᵖ
```

If ``x ∝ Uᵅ`` and ``y ∝ E²``, this function will fit the following curve:

```math
yᵅ = a₀ + a₁⋅x + a₂⋅x² + a₃⋅x³ + ⋯ + aₚ⋅xᵖ
```

The number order of the polynomial on the right hand side is given by
keyword parameter `N`.


## Examples
```@example
julia> U1 = 1.0:0.2:5.0
1.0:0.2:5.0

julia> E = sqrt.(2.0 .+ 2.0 * sqrt.(U1));

julia> U = (-1.0 .+ 0.5*E.^2).^2;

```
"""
function KingPoly(x::AbstractVector{T}, y::AbstractVector{T}, a, N) where {T}

     ya = y .^ a
   

    p1 = Polynomials.fit(x, ya, N)
    
    # Let's try to fit the polynomial

     # Fit the nonzero velocities to a polynomial of degree N

    return KingPoly{T}(p1, a)

end


function (fit::KingPoly)(x)

    ya = fit.poly(x)

    if ya < 0
        return 0 * ya
    else
        return ya ^ (1/fit.a)
    end
end

"""
 * `makekingfitfun(;a=0.45, N=4)`


This is a helper function that returns an anonymous function that is able to fit
calibration data to a modified King's Polynomia, see [`KingPoly`](@ref).


## Examples

```@example
julia> U = 1.0:0.2:5.0; E = sqrt.(2.0 .+ 2.0 * sqrt.(U));

```
"""
makekingfitfun(;a=0.45, N=4) =
    (x,y) -> KingPoly(x, y, a, N) 



    
