"""
    NDE

Author: David James, davidabraham@ucla.edu
Date: 20200829
Notes: algorithms and descriptions come from the following
    "Chapter 10: Numerical Solutions of Nonlinear Systems." Numerical Analysis,
        by Richard L. Burden et al., Cengage Learning 2016.

Contains:
- jacobian
- newtons_system
"""
module NDE

include("linalg.jl")
using .LinAlg
const la = LinAlg

export jacobian, newtons_system

"""
    jacobian(f::Function, x0::AbstractVector; ϵ::Float64=1e-10)

Calculate the jacobian of a system of equations at a `x0`.

# Arguments
- `f::Function` : the system of equations
- `x0::AbstractVector` : the point of interest for the jacobian
- `ϵ::Float64=1e-10` : a small-step taken for the derivation
"""
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

"""
    newtons_system(f::Function, x0::AbstractVector; tol::Float64=1e-5, N::Int64=50)

Approximate the solution of the nonlinear system `F`(`x`) = 0, given an initial
approximation `x`.

# Arguments
- `f::Function` : the system of equations
- `x0::AbstractVector` : the initial approximation
"""
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
