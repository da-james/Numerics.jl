"""
    BVP

Author: David James, davidabraham@ucla.edu
Date: 20200829
Notes: algorithms and descriptions come from the following
    "Chapter 11: Boundary-Value Problems for ODEs." Numerical Analysis,
        by Richard L. Burden et al., Cengage Learning 2016.

Contains:
- finite_difference
"""
module BVP

include("linalg.jl")

using .LinAlg
const la = LinAlg

export finite_difference

"""
    finite_difference(p::Function, q::Function, r::Function, a::Real, b::Real, α::Real, β::Real, N::Int64)

Approximate the solution of the boundary-value problem
    y'' = p(x)*y' + q(x)*y + r(x), for a <= x <= b, with y(a)=α and y(b)=β

# Arguments
- `p::Function` : the p(x) equation next to y'
- `q::Function` : the q(x) equation next to y
- `r::Function` : the r(x) equation
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `α::Real` : the initial conditions for the system s.t. y(a)=α
- `β::Real` : the initial conditions for the system s.t. y(b)=β
- `N:Int64` : the spacing on the interval of [a,b]
"""
function finite_difference(p::Function, q::Function, r::Function, a::Real, b::Real, α::Real, β::Real, N::Int64)

    h = (b - a) / (N + 1)

    x = a:h:b
    A = zeros(N, N + 1)

    A[1,N+1] = -h^2 * r(x[2]) + (1 + h / 2 * p(x[2])) * α
    A[N,N] = 2 + h^2 * q(x[N+1])

    for i in 2:N
        s1 = i - 1
        i1 = i + 1

        A[s1,s1+1] = -1 + h / 2 * p(x[i])
        A[s1,s1] = 2 + h^2 * q(x[i])
        A[i,s1] = -1 - h / 2 * p(x[i1])

        if i < N
            A[i,N+1] = -h^2 * r(x[i1])
        else
            A[i,N+1] = -h^2 * r(x[i1]) + (1 - h / 2 * p(x[i1])) * β
        end
    end

    w = la.crout_factorization(A)

    n = size(x)[1]
    u = zeros(n, 2)

    u[:,1] = x

    u[1,2] = α
    u[2:n-1,2] = w
    u[n,2] = β

    return u
end

end # end of module
