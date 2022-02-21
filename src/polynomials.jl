module Polynomials

export plgndr, psphrc

"""
    plgndr(x::Real, l::Int, m::Int)::Real

Computes the associated Legendre polynomial `Pᵐ_l(x)`. Here `m` and
`l` are integers satisfying `0 ≤ m ≤ l`, while `x` lies in the
range `-1 ≤ x ≤ 1`.
"""
function plgndr(x::Real, l::Int, m::Int)::Real

    P = 0

    if(m < 0 || m > l)
        println("Improper parameters, 0 ≤ m ≤ l")
        return nothing
    elseif(abs(x) > 1)
        println("Improper parameters, -1 ≤ x ≤ 1")
        return nothing
    end

    # computing Pᵐ_m
    p_mm = 1

    if(m > 0)

        x2 = sqrt((1 - x) * (1 + x))
        f = 1

        for i in 1:m
            p_mm *= -f * x2
            f += 2
        end
    end

    if(l == m)

        P = p_mm

    else

        # compute Pᵐ_m+1
        p_mm1 = x * (2*m + 1) * p_mm

        if(l == m + 1)
            P = p_mm1
        else
            # compute Pᵐ_(l > m+1)
            p_ll = 0
            for ll in m+2:l
                p_ll = (x*(2*ll - 1)*p_mm1 -(ll + m - 1)*p_mm) / (ll - m)
                p_mm = p_mm1
                p_mm1 = p_ll
            end
            P = p_ll
        end
    end

    return P

end

"""
    psphrc(θ::Real, ϕ::Real, l::Int, m::Int)::Number

Computes the associated Spherical Harmonic `Yᵐ_l(θ,ϕ)`. Here `m` and
`l` are integers satisfying `0 ≤ m ≤ l`. While (θ,ϕ) are coordinates
on the surface of a sphere.
"""
function psphrc(θ::Real, ϕ::Real, l::Int, m::Int)::Number

    fact = factorial(big(l - m)) / factorial(big(l + m))
    c = sqrt((2l + 1)/(4π) * fact)

    x = cos(θ)
    Y = c * plgndr(x, l, m) * exp(m*ϕ * im)

    if(ϕ == 0 || m == 0)
        return real(Y)
    else
        return Y
    end
end

end # module
