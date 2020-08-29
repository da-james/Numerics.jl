module BVP

include("linalg.jl")

using .LinAlg

const la = LinAlg

function linear_difference(p::Function, q::Function, r::Function, a::Real, b::Real, α::Real, β::Real, N::Int64)

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
