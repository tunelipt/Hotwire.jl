using Polynomials


struct HWKingFit
    A::Float64
    B::Float64
    n::Float64
end
Base.broadcastable(fit::HWKingFit) = Ref(fit)

function HWKingFit(E::AbstractVector, U::AbstractVector, ϵ=1e-8, maxiter=200)

    A, B, n = king_fit(E, U, ϵ, maxiter)

    return HWKingFit(A, B, n)
end

function (fit::HWKingFit)(E)

    E² = E*E
    if E² < fit.A
        return zero(E)
    else
        return (  (E²-fit.A) / fit.B )^(1/fit.n)
    end
end

    
struct KingPoly
    poly::Polynomial{Float64, :E²}
    n::Float64
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

\$Uⁿ = a₀ + a₁⋅E² + a₂⋅(E²)² + a₃⋅(E²)³ + ⋯ + aₚ⋅(E²)ᵖ\$

If `n` is provided, a linear fit is carried out. If it is not provided,
an initial guess is used (`n = 0.4`) a Newton-like iteration is used to improve
the values.

The number of terms on the right hand side is given by keyword parameter `N`.

"""
function KingPoly(E::AbstractVector, U::AbstractVector, n; N=5)

    Ur = U .^ n
    E2 = E .^ 2

    # Let's try to fit the polynomial

     # Fit the nonzero velocities to a polynomial of degree N
    coeffs = poly_fit(E2, Ur, N) 

    return KingPoly(Polynomial(coeffs, :E²), n)

end

function KingPoly(E::AbstractVector, U::AbstractVector; N=5, err=1e-9, maxiter=200)
    n = 0.4
    Ur = U .^ n
    E2 = E .^ 2

    # Let's try to fit the polynomial

     # Fit the nonzero velocities to a polynomial of degree N
    coeffs = poly_fit(E2, Ur, N) 
    a = [coeffs; n]
    #return a
    a1, conv, iter = nonlinear_fit([E U], king_poly_aux, a, err, maxiter)

    poly = Polynomial(a1[1:end-1], :E²)

    return KingPoly(poly, a1[end])
    
end

function king_poly_aux(x, a)
    E2 = x[1]^2
    U = x[2]
    Ur = evalpoly(E2, a) - a[end]*E2^(length(a)-1)
    return U^a[end] - Ur
end


function (fit::KingPoly)(E)

    U1 = fit.poly(E*E)

    if U1 < 0
        return zero(U1)
    else
        return U1 ^ (1/fit.n)
    end
end

