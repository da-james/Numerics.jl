module Polynomials

export plgndr

"""
    plgndr(x::Real, l::Int, m::Int)

Computes the associated Legendre polynomail `Pᵐ_l(x)`. Here `m` and
`l` are integers satisfying `0 ≤ m ≤ l`, while `x` lies in the
range `-1 ≤ x ≤ 1`.
"""
function plgndr(x::Real, l::Int, m::Int)

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

end # module
