module BVP

include("linalg.jl")

using .LinAlg

const la = LinAlg

function linear_difference(p::Function, q::Function, r::Function, a::Real, b::Real, α::Real, β::Real, N::Int64)

    h = (b - a) / N

    x = a:h:b
    println(size(x))
    println(x[1])
    println(x[2])
    println(x[end])
    A = zeros(N, N + 1)

    A[1,N+1] = -h^2 * r(x[1]) + (1 + h / 2 * p(x[1])) * α
    A[N,N] = 2 + h^2 * q(x[N])

    for i in 2:N
        s1 = i-1
        A[s1,s1+1] = -1 + h / 2 * p(x[s1])
        A[i-1,i-1] = 2 + h^2 * q(x[i-1])
        A[i,i-1] = -1 - h / 2 * p(x[i-1])

        if i < N
            A[1,i] = -h^2 * r(x[i])
        else
            A[1,i] = -h^2 * r(x[i]) + (1 - h / 2 * p(x[i])) * β
        end
    end

    w = la.crout_factorization(A)
    println(size(w))

    return (x, w)
end

end # end of module
