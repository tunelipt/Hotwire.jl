struct HWKingFit
    A::Float64
    B::Float64
    n::Float64
end
Base.broadcastable(fit::HWKingFit) = Ref(cal)

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

    
