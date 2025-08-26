using Polynomials

"""
 * `linearfit(y::AbstractVector{T}, ϕ::AbstractVector{T})`
 * `linearfit(y::AbstractVector{T}, ϕ₁::AbstractVector{T}, , ϕ₁::AbstractVector{T})`
 * `linearfit(y::AbstractVector{T}, ϕ...)`
 * `linearfit(y::AbstractVector{T}, ϕ::AbstractMatrix{T})`

This function uses the the least squares principle to fit a linear model to data.

The first form fits

`y[i] = a ⋅ ϕ[i]`

returning coefficient `a`

The second form fits

`y[i] = a₁ ⋅ ϕ₁[i] + a₂ ⋅ ϕ₂[i]`

returning the tuple `(a₁, a₂)`.

The third form fits

`y[i] = a₁ ⋅ ϕ₁[i] + a₂ ⋅ ϕ₂[i] + ⋯ + aₖ ⋅ ϕₖ[i]`

returning the vector `[a₁, a₂, …, aₖ]`.

Finally the last form where `ϕ` is a matrix, the function fits

`y[i] = a₁ ⋅ ϕ[i,1] + a₂ ⋅ ϕ[i,2] + ⋯ + aₖ ⋅ ϕ[i,k]`

returning the vector `[a₁, a₂, …, aₖ]`.

 

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

function linearfit(y::AbstractVector{T}, ϕ...) where {T}
    Nϕ = length(ϕ)
    N = length(y)

    for i in 1:Nϕ
        @assert length(ϕ[i]) == N "ϕ $i has wrong length!"
    end

    A = zeros(T, Nϕ, Nϕ)
    b = zeros(T, Nϕ)
    for k in 1:Nϕ
        ϕₖ = ϕ[k]
        for j in k:Nϕ
            ϕⱼ = ϕ[j]
            AA = zero(T)
            for (ϕₖⁱ,ϕⱼⁱ) in zip(ϕₖ, ϕⱼ)
                AA += ϕₖⁱ * ϕⱼⁱ
            end
            A[k,j] = AA
            A[j,k] = AA
        end
        bb = zero(T)
        for (ϕₖⁱ,yᵢ) in zip(ϕₖ,y)
            bb += yᵢ * ϕₖⁱ
        end
        b[k] = bb
    end

    return A\b
end


function linearfit(y::AbstractVector{T}, ϕ::AbstractMatrix{T}) where {T}

    @assert length(y) == size(ϕ,1)

    Nϕ = size(ϕ,2)
    N = length(y)
    
    A = zeros(T,Nϕ,Nϕ)
    b = zeros(T,Nϕ)

    for k in 1:Nϕ
        for j in k:Nϕ
            AA = zero(T)
            for i in 1:N
                AA += ϕ[i,k] * ϕ[i,j]
            end
            A[k,j] = AA
            A[j,k] = AA
        end
        bb = zero(T)
        for i in 1:N
            bb += y[i]*ϕ[i,k]
        end
        b[k] = bb
    end

    return A\b
end

function polyfit1(y::AbstractVector{T}, x::AbstractVector{T}, P::Integer) where {T}
    N = length(y)
    @assert N == length(x)
    A = zeros(T, P+1, P+1)
    b = zeros(T, P+1)

    for k in 1:P+1
        k₁ = k-1
        for j in k:P+1
            j₁ = j-1
            a = zero(T)
            for i in 1:N
                a += x[i]^k₁ * x[i]^j₁
            end
            A[k,j] = a
            A[j,k] = a
        end
        bb = zero(T)
        for i in 1:N
            bb += y[i] * x[i]^k₁
        end
        b[k] = bb
    end
    return A\b
end

struct KingLaw{T}
    A::T
    B::T
    n::T
end
Base.broadcastable(fit::KingLaw) = Ref(fit)

"""
 * `KingLaw(A,B,n)`
 * `KingLaw(E,U,n)`
 * `KingLaw(E,U, n0=0.4)`

Fits hotwire anemometry calibration data to modified King's law. It finds the 
coefficients `A`, `B` and `n` in equation

```math
E² = A + B ⋅ Uⁿ
```

If `n` is provided as a positional argument, the value provided for `n` is fixed.
If the argument `n` is not provided or the keyword argument `n0`, `n` will be fitted using nonlineas least squares with the keyword argument `n0` used as initial value.

If there are problems in convergence, keyword parameter `r` is a subrelaxation
factor and can help in convergence. The convergence is checked when the step in
power law exponent `n` is less than the maximum convergence error `err`.
The maximum number of iterations is given by keyword parameter `maxiter`.
If convergence fails, an error is thrown.


A more general and often more accurate fit function is implemented
in [`KingPoly`](@ref).


## Examples
```@example
ulia> U = 1.0:0.2:5.0
1.0:0.2:5.0

julia> E = sqrt.(2.0 .+ 2.0 * sqrt.(U));

julia> KingLaw(E,U,0.5)
KingLaw{Float64}(2.0000000000000258, 1.999999999999985, 0.5)

julia> KingLaw(E,U, n0=0.4)
KingLaw{Float64}(2.0000000000000053, 1.9999999999999951, 0.5000000000000009)
```

"""
function KingLaw(E::AbstractVector{T}, U::AbstractVector{T}, n) where {T}

    E² = E .* E
    Uⁿ = U .^ n

    A,B = linearfit(E², fill(one(T),length(E)), Uⁿ)
    return KingLaw(A,B,n)

end

function KingLaw(E::AbstractVector{T}, U::AbstractVector{T};
                 n=one(T)/2, fixed=false, err=1e-8, maxiter=200, r=1.0) where {T}
    if fixed
        return KingLaw(E, U, n)
    else
        n0 = n
    end
    
    M = length(E)
    @assert length(U) == M  "E and U should have the same length"

    E² = E .* E
    lnU = log.(U)
    y = zeros(T, M)

    o = ones(T,M)
    ∂a₀ = o
    ∂a₁ = E²
    ∂n  = zeros(T,M)
    n = n0
    Uⁿ = U .^ n
    a₀,a₁ = linearfit(Uⁿ, o, E²)
    ϵ = zeros(T,M)

    
    for iter in 1:maxiter
        Uⁿ .= U .^ n
        y .=  a₀ .+ a₁ .* E²
        ϵ .=  (y .- Uⁿ) 
        ∂n .= Uⁿ .* lnU

        δ = linearfit(ϵ, ∂n, ∂a₀, ∂a₁)

        if abs(δ[1]) < err
            n  +=  δ[1]
            a₀ -=  δ[2]
            a₁ -=  δ[3]
            
            return KingLaw(-a₀/a₁, 1/a₁, n)
        else
            n  +=  r*δ[1]
            a₀ -=  r*δ[2]
            a₁ -=  r*δ[3]
        end
    end

    # If we got here, the algorithm failed to converge
    error("`KingLaw` failed to converge to residual $err in $maxiter iterations.")
    
end

errfun(E,U,fit) = E*E - (fit.A + fit.B * U^fit.n)


function (fit::KingLaw{T})(E) where {T}

    E² = E*E
    if E² < fit.A
        return zero(T)
    else
        return (  (E²-fit.A) / fit.B )^(1/fit.n)
    end
end

    
struct KingPoly{T}
    poly::Polynomial{T, :E²}
    n::T
end

Base.broadcastable(fit::KingPoly) = Ref(fit)

"""
`KingPoly(E, U, n; N=5)`

Fits a King Polynomial to the calibration curve of a thermal anemometer.

King's law has strong connection to the physics of the problem bu it
is often not accurate enough. On the other hand, a simple polynomial
between `E` and `U` does work well but has all the problems associated
with polynomials. A possibility is to get the best of both worlds:
expand King's law!

The proposed approximation is:

```math
Uⁿ = a₀ + a₁⋅E² + a₂⋅(E²)² + a₃⋅(E²)³ + ⋯ + aₚ⋅(E²)ᵖ
```
The number order of the polynomial on the right hand side is given by
keyword parameter `N`.

If the `n` positional argument is provided, a linear fit is carried out.
If it is not provided, a nonlinear least squares Newton iteration will be
 used with the keyword argument `n0` used as initial guess.

If there are problems in convergence, keyword parameter `r` is a subrelaxation
factor and can help in convergence. The convergence is checked when the step in
power law exponent `n` is less than the maximum convergence error `err`.
The maximum number of iterations is given by keyword parameter `maxiter`.
If convergence fails, an error is thrown.


The function [`makekingpolyfitfun`](@ref) returns an anonymous function that has two
arguments `E` and `U` that can be used to fit the data.

See also [`KingLaw`](@ref) for the conventional modified King's law fitting function.
## Examples
```@example
julia> U1 = 1.0:0.2:5.0
1.0:0.2:5.0

julia> E = sqrt.(2.0 .+ 2.0 * sqrt.(U1));

julia> U = (-1.0 .+ 0.5*E.^2).^2;

julia> KingPoly(E,U, n=0.5,N=1)
KingPoly{Float64}(Polynomial(-0.9999999999999994 + 0.5*E²), 0.5)

julia> KingPoly(E,U,N=1,n0=0.4)
KingPoly{Float64}(Polynomial(-0.9999999999999992 + 0.4999999999999998*E²), 0.49999999999999983)

julia> KingPoly(E,U,N=3,n0=0.4)
KingPoly{Float64}(Polynomial(-1.0000000000000018 + 0.49999999999998945*E² + 3.1445640847865035e-15*E²^2 - 9.783317428355699e-17*E²^3), 0.5000000000000097)

```
"""
function KingPoly(E::AbstractVector{T}, U::AbstractVector{T}, n, N) where {T}

    Ur = U .^ n
    E2 = E .^ 2

    p1 = Polynomials.fit(E2, Ur, N, var=:E²)
    
    # Let's try to fit the polynomial

     # Fit the nonzero velocities to a polynomial of degree N

    return KingPoly{T}(p1, n)

end

function KingPoly(E::AbstractVector{T}, U::AbstractVector{T};
                  N=5, n=one(T)/2, fixed=true, err=sqrt(eps(T)),
                  maxiter=200, r=one(T)) where {T}
    if fixed
        return KingPoly(E, U, n, N)
    else
        n0 = n
    end
    
    M = length(E)
    @assert length(U) == M  "E and U should have the same length"

    E² = E .* E
    lnU = log.(U)
    y = zeros(T, M)

    o = ones(T,M)
    ∂ = zeros(T,M,N+2)
    for k in 0:N
        for i in 1:M
            ∂[i,k+2] = E²[i] .^ k
        end
    end
    
    ∂n  = zeros(T,M)
    n = n0
    Uⁿ = U .^ n
    p = Polynomials.fit(E², Uⁿ, N, var=:E²)
    a = coeffs(p)
    
    ϵ = zeros(T,M)

    
    for iter in 1:maxiter
        Uⁿ .= U .^ n
        y .=  p.(E²)
        ϵ .=  (y .- Uⁿ) 
        ∂n .= Uⁿ .* lnU
        for i in 1:M
            ∂[i,1] = ∂n[i]
        end
        
        δ = linearfit(ϵ, ∂)

        if abs(δ[1]) < err
            n  +=  δ[1]
            for k in 2:N+2
                a[k-1] -=  δ[k]
            end
            return KingPoly(p, n)
        else
            n  +=  r*δ[1]
            for k in 2:N+2
                a[k-1] -=  r*δ[k]
            end
        end
    end
    # If we got here, the algorithm failed to converge
    error("`KingPoly` failed to converge to residual $err in $maxiter iterations.")

end

function (fit::KingPoly{T})(E) where {T}

    U1 = fit.poly(E*E)

    if U1 < 0
        return zero(T)
    else
        return U1 ^ (1/fit.n)
    end
end

"""
 * `make_king_fitfun()`
 * `make_king_fitfun(0.5)`

This is a helper function that returns an anonymous function that is able to fit
calibration data to a modified King's Law, see [`KingLaw`](@ref).

When used without positional arguments, it uses nonlinear least squares iteration.
If positional argument `n` is provided, the power law exponent `n` is fixed.

## Examples

```@example
julia> U = 1.0:0.2:5.0; E = sqrt.(2.0 .+ 2.0 * sqrt.(U));

julia> f = make_king_fitfun()
#34 (generic function with 1 method)

julia> f(E,U)
KingLaw{Float64}(2.0000000000000018, 1.9999999999999982, 0.5000000000000003)
```
"""
makekingfitfun(;n0=0.45, err=1e-8, maxiter=200, r=1.0) =
    (E,U) -> KingLaw(E,U, n0=n0, err=err, maxiter=maxiter, r=r)
makekingfitfun(n) = (E,U) -> KingLaw(E,U,n)



struct PowerPoly{T}
    poly::Polynomial{T, :E²}
    n::T
end
Base.broadcastable(fit::PowerPoly) = Ref(fit)

function PowerPoly(E²::AbstractVector{T}, U::AbstractVector{T}; n=0.5, N=5) where {T<:AbstractFloat}
    Uⁿ = U .^ n
    p1 = Polynomials.fit(E², Uⁿ, N, var=:E²)
    
    return PowerPoly(p1, n)
end


function (fit::PowerPoly{T})(E²) where {T}

    U1 = fit.poly(E²)

    if U1 < 0
        return zero(T)
    else
        return U1 ^ (1/fit.n)
    end
end

    
