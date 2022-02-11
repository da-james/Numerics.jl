"""
    ODE

Author: David James, davidabraham@ucla.edu
Date: 20200829
Notes: algorithms and descriptions come from the following
    "Chapter 5: Initial-Value Problems for ODEs." Numerical Analysis,
        by Richard L. Burden et al., Cengage Learning 2016.

Contains:
- eulers
- rk4
- trapezoid
"""
module ODE

using ..Numerics: norm

export eulers, rk4, trapezoid

"""
    eulers(f::Function, α::Real, a::Real, b::Real, N::Int64)

Approximate the solution of the initial-value problem (IVP)
    y' = f(t,y), a <= t <= b, y(a) = α
at (N + 1) equally spaced steps in the interval [a,b].
NOTE: eulers is meant is a education tool for a more robust
      solver, refer to the `rk4` method.

# Arguments
- `f::Function` : the ODE
- `α::Real` : the initial condition such that y(a) = α
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `N:Int64` : the number of steps on [a,b]
"""
function eulers(f::Function, α::Real, a::Real, b::Real, N::Int64)

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
    rk4(f::Function, α::AbstractArray, a::Real, b::Real, N::Int64, p)

Approximate the solution of the mth-order system of a first-order IVP

    u'ⱼ = fⱼ(t, u₁, u₂, …, uₘ, p), a <= t <= b, with uⱼ(a) = αⱼ

for `j = 1, 2, …, m` at `(N + 1)` equally spaced numbers in the
interval `[a,b]`. The `p` argument are any parameter values needed for
the function.

# Arguments
- `f::Function`      : the system of equations of ODEs
- `α::AbstractArray` : the initial conditions for the system
- `a::Real`          : the left-sided endpoint
- `b::Real`          : the right-sided endpoint
- `N:Int64`          : the number of steps on [a,b]
- `p`                : parameters for `f(t, u, p)`
"""
function rk4(f::Function, α::AbstractArray, a::Real, b::Real,
             N::Int64, p)

    n1 = N + 1
    m = size(α)[1]

    u = zeros(n1, m + 1)

    h = (b - a) / N
    u[1,1] = a
    u[1,2:end] = α

    for i in 2:n1
        t = u[i-1,1]
        v = u[i-1,2:end]

        k1 = h * f(t, v, p)
        k2 = h * f(t + h / 2, v .+ k1 / 2, p)
        k3 = h * f(t + h / 2, v .+ k2 / 2, p)
        k4 = h * f(t + h, v .+ k3, p)

        u[i,2:end] = v .+ (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) ./ 6
        u[i,1] = a + (i - 1) * h
    end

    return u
end

"""
    rkf45(f::Function, α::AbstractArray, a::Real, b::Real,
          h::Real, h_min::Real, p; tol::Real=1e-6)

Approximate the solution of the mth-order system of a first-order IVP

    u'ⱼ = fⱼ(t, u₁, u₂, …, uₘ, p), a <= t <= b, with uⱼ(a) = αⱼ

for `j = 1, 2, …, m` at `(N + 1)` equally spaced numbers in the
interval `[a,b]`. The `p` argument are any parameter values needed for
the function.

# Arguments
- `f::Function`      : the system of equations of ODEs
- `α::AbstractArray` : the initial conditions for the system
- `a::Real`          : the left-sided endpoint
- `b::Real`          : the right-sided endpoint
- `N:Int64`          : the number of steps on [a,b]
- `p`                : parameters for `f(t, u, p)`
"""
function rkf45(f::Function, α::AbstractArray, a::Real, b::Real,
               h::Real, h_min::Real, p; tol::Real=1e-6)

    l = length(α)
    u = zeros(1, l + 1)

    u[1,1] = a
    u[1,2:end] = α

    t_k = a
    y_k = α

    while t_k <= b

        k1 = h*f(t_k, y_k)
        k2 = h*f(t_k + 1/4*h  , y_k .+ 1/4*k1)
        k3 = h*f(t_k + 3/8*h  , y_k .+ 3/32*k1      .+ 9/32*k2)
        k4 = h*f(t_k + 12/13*h, y_k .+ 1932/2197*k1 .- 7200/2197*k2 .+ 7296/2197*k3)
        k5 = h*f(t_k + h      , y_k .+ 439/216*k1   .- 8*k2         .+ 3680/513*k3  .- 845/4104*k4)
        k6 = h*f(t_k + 1/2*h  , y_k .- 8/27*k1      .+ 2*k2         .- 3544/2565*k3 .+ 1859/4104*k4 .- 11/40*k5)

        y_k1 = y_k .+ 25/216*k1 .+ 1408/2565*k3  .+ 2197/4101*k4   .- 1/5*k5
        z_k1 = y_k .+ 16/135*k1 .+ 6656/12825*k3 .+ 28561/56430*k4 .- 9/50*k5 .+ 2/55*k6

        R = 1/4 * norm(z_k1 .- y_k1)

        if(R <= tol)
            t_k += h
            y_k = y_k1
            v_k = vcat([t], y_k)
            u = vcat(u, v_k)
        else
            δ = (tol / (2*R))^(1/4)
            h *= δ
            if(h < h_min) h = h_min end
        end
    end

    return u

end


"""
    trapezoid(f::Function, df::Function, α::Real, a::Real, b::Real,
              N::Int64; tol::Float64=1e-5, M::Int64=10)

Approximate the solution of the stiff initial-value problem (IVP)
    y' = f(t,y), a <= t <= b, y(a) = α
at (N + 1) equally spaced steps in the interval [a,b].
NOTE: The trapezoid isn't fully implemented to take in parameters
      or handle system of equations. Refer to the method `rk4` for
      those needs, but it may be unstable for stiff equations.

# Arguments
- `f::Function` : the ODE as f(t,y)
- `df::Function` : the derivative of f(t,y)
- `α::Real` : the initial condition such that y(a) = α
- `a::Real` : the left-sided endpoint
- `b::Real` : the right-sided endpoint
- `N:Int64` : the number of steps on [a,b]
"""
function trapezoid(f::Function, df::Function, α::Real, a::Real, b::Real,
                   N::Int64; tol::Float64=1e-5, M::Int64=10)

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
