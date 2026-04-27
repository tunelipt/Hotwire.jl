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



    




"""
 * `linearfit(y::AbstractVector{T}, ϕ::AbstractVector{T})`
 * `linearfit(y::AbstractVector{T}, ϕ₁::AbstractVector{T}, , ϕ₁::AbstractVector{T})`

This function uses the the least squares principle to fit a linear model to data.

The first form fits

`y[i] = a ⋅ ϕ[i]`

returning coefficient `a`

The second form fits

`y[i] = a₁ ⋅ ϕ₁[i] + a₂ ⋅ ϕ₂[i]`

returning the tuple `(a₁, a₂)`.


 

"""
function linearfit(y::AbstractVector{T}, ϕ::AbstractVector{T}) where {T}
    @assert length(y) == length(ϕ) "Incompatible sizes!"
    
    A = zero(T)
    b = zero(T)
    for (yᵢ, ϕᵢ) in zip(y,ϕ)
        A += ϕᵢ*ϕᵢ
        b += yᵢ*ϕᵢ
    end
    return b / A
end

function linearfit(y::AbstractVector{T},
                   ϕ₁::AbstractVector{T}, ϕ₂::AbstractVector{T}) where {T}

    @assert length(y) == length(ϕ₁) == length(ϕ₂) "Incompatible sizes!"
    
    A₁₁ = zero(T)
    A₁₂ = zero(T)
    A₂₂ = zero(T)
    b₁ = zero(T)
    b₂ = zero(T)

    for i in eachindex(y)
        yᵢ = y[i]
        ϕ₁ᵢ = ϕ₁[i]
        ϕ₂ᵢ = ϕ₂[i]
                
        b₁ += ϕ₁ᵢ * yᵢ
        b₂ += ϕ₂ᵢ * yᵢ

        A₁₁ += ϕ₁ᵢ * ϕ₁ᵢ
        A₂₂ += ϕ₂ᵢ * ϕ₂ᵢ
        A₁₂ += ϕ₁ᵢ * ϕ₂ᵢ
    end

    a₂ = (b₂ - A₁₂/A₁₁ * b₁) / (A₂₂ - A₁₂*A₁₂/A₁₁)
    a₁ = (b₁  - A₁₂ * a₂) /  A₁₁

    return a₁, a₂
        
end
