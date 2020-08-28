module Nonlinear

include("linalg.jl")

using .LinAlg

const la = LinAlg

function jacobian(f::Function, x0::AbstractVector; ϵ::Float64=1e-10)

    n = size(x0)[1]
    J = zeros(n, n)

    fx = f(x0)
    x_perb = copy(x0)

    for i in 1:n
        x_perb[i] += ϵ
        J[:,i] = (f(x_perb) - fx) / ϵ
        x_perb[i] = x0[i]
    end

    return J
end

function newtons_system(f::Function, x0::AbstractVector; tol::Float64=1e-5, N::Int64=50)
    k = 1

    while k <= N
        J = jacobian(f, x0)
        fx = -1 .* f(x0)

        y = la.gauss_sidel_method(J, fx, x0)
        x0 += y

        if la.norm(y) < tol
            return x0
        end

        k += 1
    end

    println("Maximum number of iterations exceeded")
    return nothing
end

end # end of module
