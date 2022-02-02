"""
    ODE

Author: David James, davidabraham@ucla.edu
Date: 20200829
Notes: algorithms and descriptions come from the following
    "Chapter 5: Initial-Value Problems for ODEs." Numerical Analysis,
        by Richard L. Burden et al., Cengage Learning 2016.

Contains:
- eulers_method
- rk4_method
- rk4_system
- trapezoid_method
"""
module ODE

export rk4_method, rk4_system, trapezoid_method

"""
    eulers_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

Approximate the solution of the initial-value problem (IVP)
    y' = f(t,y), a <= t <= b, y(a) = α
at (N + 1) equally spaced numbers in the interval [a,b].

# Arguments
- `f::Function` : the ODE
- `α::Real` : the initial condition such that y(a) = α
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `N:Int64` : the spacing on the interval of [a,b]
"""
function eulers_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

    n1 = N + 1
    u = zeros(n1,2)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        u[i,2] = u[i-1,2] + h * f(u[i-1,1], u[i-1,2])
        u[i,1] = a + (i - 1) * h
    end

    return u
end

"""
    rk4_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

Approximate the solution of the initial-value problem (IVP)
    y' = f(t,y), a <= t <= b, y(a) = α
at (N + 1) equally spaced numbers in the interval [a,b].

# Arguments
- `f::Function` : the ODE
- `α::Real` : the initial condition such that y(a) = α
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `N:Int64` : the spacing on the interval of [a,b]
"""
function rk4_method(f::Function, α::Real, a::Real, b::Real, p, N::Int64)

    n1 = N + 1
    u = zeros(n1,2)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        t = u[i-1,1]
        w = u[i-1,2]

        k1 = h * f(t, w, p)
        k2 = h * f(t + h / 2, w + k1 / 2, p)
        k3 = h * f(t + h / 2, w + k2 / 2, p)
        k4 = h * f(t + h, w + k3, p)

        u[i,2] = w + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        u[i,1] = a + (i - 1) * h
    end

    return u
end

"""
    rk4_system(f::Function, α::AbstractArray, a::Real, b::Real, N::Int64)

Approximate the solution of the mth-order system of first-order IVP
    u'_j = f_j(t, u1, u2, ..., u_m), a <= t <= b, with u_j(a) = α_j
for j = 1, 2, ..., m at (N + 1) equally spaced numbers in the interval [a,b].

# Arguments
- `f::Function` : the system of equations of ODEs
- `α::Real` : the initial conditions for the system
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `N:Int64` : the spacing on the interval of [a,b]
"""
function rk4_system(f::Function, α::AbstractArray, a::Real, b::Real, p, N::Int64)

    n1 = N + 1
    m = size(α)[1]

    u = zeros(n1, m + 1)

    h = (b - a) / N
    u[1,1] = a
    u[1,2:end] = α

    for i in 2:n1
        t = u[i-1,1]

        v = u[i-1,2:end]

        k1 = h * f(t, v)
        k2 = h * f(t + h / 2, v .+ k1 / 2)
        k3 = h * f(t + h / 2, v .+ k2 / 2)
        k4 = h * f(t + h, v .+ k3)

        u[i,2:end] = v .+ (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) ./ 6
        u[i,1] = a + (i - 1) * h
    end

    return u
end

"""
    trapezoid_method(f::Function, df::Function, α::Real, a::Real, b::Real, N::Int64; tol::Float64=1e-5, M::Int64=10)

Approximate the solution of the stiff initial-value problem (IVP)
    y' = f(t,y), a <= t <= b, y(a) = α
at (N + 1) equally spaced numbers in the interval [a,b].

# Arguments
- `f::Function` : the ODE as f(t,y)
- `df::Function` : the derivative of f(t,y)
- `α::Real` : the initial condition such that y(a) = α
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `N:Int64` : the spacing on the interval of [a,b]
"""
function trapezoid_method(f::Function, df::Function, α::Real, a::Real, b::Real, N::Int64; tol::Float64=1e-5, M::Int64=10)

    n1 = N + 1
    u = zeros(n1,2)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        ti = u[i-1,1]
        wi = u[i-1,2]

        k1 = wi + h / 2 * f(ti, wi)
        w0 = k1
        j = 1
        FLAG = 0

        while FLAG == 0
            num = w0 - (h / 2 * f(ti + h, w0)) - k1
            dem = 1 - h / 2 * df(ti + h, w0)

            w = w0 - num / dem

            if abs(w - w0) < tol
                FLAG = 1
            else
                j += 1
                w0 = w
                if j > M
                    println("The maximum number of iterations exceeded")
                    return nothing
                end
            end
            u[i, 2] = w
        end

        u[i, 1] = a + (i - 1) * h
    end

    return u
end

end # end of module
