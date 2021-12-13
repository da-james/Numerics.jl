module Methods

function bisection_method(f::Function, a::Real, b::Real; tol::Float64=1e-5, N::Int64=50)

    # first iteration
    n = 1
    # left bound
    fa = f(a)

    a_vals = []
    p_vals = []
    b_vals = []

    # while less than max iterations
    while n <= N
        # mid point of bounds
        mid = (b - a) / 2
        p = a + mid
        fp = f(p)

        # updating arrays
        push!(a_vals, a)
        push!(p_vals, p)
        push!(b_vals, b)

        # if mid point is less than tolerance
        # then return p
        if fp == 0 || abs(mid) < tol
            println("iterations: ", n)
            return (a_vals, p_vals, b_vals)
        end

        # increment iteration
        n += 1

        # if fa*fp is positive
        # then set left bound to p
        if fa * fp > 0
            a = p
            fa = fp
        else
            # else set right bound to p
            b = p
        end
    end

    println("Method failed after " * string(N) * "iterations.")
end

function newtons_method(f::Function, g::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)
    n = 1

    p_vals = []

    # while less than max iterations
    while n <= N
        # computing p_n
        p = p0 - f(p0) / g(p0)

        push!(p_vals, p)
        # check difference between p and p0
        if abs(p - p0) < tol
            return p_vals
        end

        # iterate
        n += 1
        # adjust p0 to new p
        p0 = p
    end

    println("Method failed after " * string(N) * "iterations.")
end

function secant_method(f::Function, p0::Real, p1::Real; tol::Float64=1e-5, N::Int64=50)
    n = 2
    q0 = f(p0)
    q1 = f(p1)

    p_vals = []

    while n <= N
        p = p1 - (q * (p1 - p0)) / (q1 - q0)

        push!(p_vals, p)
        if abs(p - p1) < tol
            return p_vals
        end

        n += 1

        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
    end

    println("Method failed after " * string(N) * "iterations.")
end

function nevilles_method(x0::Real, x::AbstractVector, q::AbstractVector)

    n = size(x)[1]

    for i in n:-1:1
        for j in 1:i
            num = (x0 - x[j]) * q[j+1] - (x0 - x[j+n+1-i]) * q[j]
            dem = x[j+n+1-i] - x[j]
            q[j] = num / dem
        end
    end

    return q
end

function newtons_difference(x::AbstractVector, q::AbstractVector)

    n = size(x)[1]
    f = zeros(n, n)
    f[:,1] = q

    for i in 2:n
        for j in i:n
            f[j,i] = (f[j,i-1] - f[j-1,i-1]) / (x[j] - x[j-i+1])
        end
    end

    return f
end

function cubic_spline(x::AbstractVector, a::AbstractVector)

    n = size(x)[1]

    h = zeros(n-1)
    α = zeros(n-1)

    for i in 1:n-1
        h[i] = x[i+1] - x[i]
        if i >= 2
            α[i] = (3 / h[i]) * (a[i+1] - a[i]) - (3 / h[i-1]) * (a[i] - a[i-1])
        end
    end

    # constructing tridiagonal linear system
    l = zeros(n)
    μ = zeros(n)
    z = zeros(n)

    b = zeros(n)
    c = zeros(n)
    d = zeros(n)

    l[1] = 1

    # solving out tridiagonal system
    # similar to crout_factorization in LinAlg module
    for i in 2:n-1
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end

    l[n] = 1

    for j in n-1:-1:1
        c[j] = z[j] - μ[j] * c[j+1]
        b[j] = (a[j+1] - a[j]) / h[j] - (c[j+1] + 2 * c[j]) * (h[j] / 3)
        d[j] = (c[j+1] - c[j]) / (3 * h[j])
    end

    return a, b, c, d
end

function eulers_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

    n1 = N + 1
    u = zeros(n1,3)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        u[i,3] = h * f(u[i-1,1], u[i-1,2])
        u[i,2] = u[i-1,2] + u[i,3]
        u[i,1] = a + (i - 1) * h
    end

    return u
end

function rk4_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

    n1 = N + 1
    u = zeros(n1,6)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        t = u[i-1,1]
        w = u[i-1,2]

        u[i,3] = h * f(t, w)
        u[i,4] = h * f(t + h / 2, w + u[i,3] / 2)
        u[i,5] = h * f(t + h / 2, w + u[i,4] / 2)
        u[i,6] = h * f(t + h, w + u[i,5])

        u[i,2] = w + (u[i,3] + 2 * u[i,4] + 2 * u[i,5] + u[i,6]) / 6
        u[i,1] = a + (i - 1) * h
    end

    return u
end

function fft(x::AbstractVector)

    N = size(x)[1]

    if N > 2
        x_odd = fft(x[1:2:N])
        x_even = fft(x[2:2:N])
    else
        x_odd = x[1]
        x_even = x[2]
    end

    n = 0:N-1
    half = div(N,2)
    factor = exp.(-2im * pi * n / N)

    return vcat(x_odd .+ x_even .* factor[1:half],
                x_odd .- x_even .* factor[1:half])
end

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
