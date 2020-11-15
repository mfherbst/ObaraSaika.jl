using SpecialFunctions


@doc raw"""
Evaluates the Boys function ``Fₘ(x) = ∫₀¹ t²ᵐ \exp(-xt²) dt``.
`IND==0` gives around 14 digits, `IND==1` around 6 and `IND==3` around 3.
"""
function boys(m::Int, x::T, IND=0) where T <: Real
    x < 0 && throw(DomainError(x, "`x` must be nonnegative."))
    if x > 0
        # gamma_inc returns tuple of incomplete Γ function ratios,
        # we only care about the first returned result
        Px, _ = gamma_inc(m + 1 / T(2), x, IND)
        1 / 2x^(m + 1 / T(2)) * gamma(m + 1 / T(2)) * Px
    else
        1 / T(2m + 1)
    end
end


function doublefactorial(n::Integer)
    n < 0 && throw(DomainError(n, "`n` must be nonnegative."))
    start = iseven(n) ? 2 : 1
    f::typeof(n*n) = 1
    for i::typeof(n*n) = start:2:n
        f *= i
    end
    return f
end
