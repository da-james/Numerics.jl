"""
    Root

Author: David James, davidabraham@ucla.edu
Date: 20200824
Notes: algorithms and descriptions come from the following
    "Chapter 2: Solutions of Equations in One Variable." Numerical Analysis, by
        Richard L. Burden et al., Cengage Learning 2016.

Contains:
- bisection_method
- fixed_point
- newtons_method
- secant_method
"""
module Root

export bisection_method, fixed_point, newtons_method, secant_method

"""
    bisection_method(f::Function, a::Real, b::Real; tol::Float64=1e-5, N::Int64=50)

Find a solution to f(x) = 0 given the continuous function `f` on the interval
[`a`,`b`], where f(`a`) and f(`b`) have opposite signs.

# Arguments
- `f::Function` : the continuous function given
- `a::Real` : the start of the interval
- `b::Real` : the end of the interval
- `tol::Float64=1e-5` : the tolerance for error when finding x such that f(x) = 0
- `N::Int64=50` : the max number of iterations to perform
"""
function bisection_method(f::Function, a::Real, b::Real; tol::Float64=1e-5, N::Int64=50)

    # first iteration
    n = 1
    # left bound
    fa = f(a)

    # while less than max iterations
    while n <= N
        # mid point of bounds
        mid = (b - a) / 2
        p = a + mid
        fp = f(p)

        # if mid point is less than tolerance
        # then return p
        if fp == 0 || abs(mid) < tol
            return p
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
    return nothing
end

"""
    fixed_point(f::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)

Find a solution p = f(p) given an initial approximation `p0`.

# Arguments
- `f::Function` : the continuous function given
- `p0::Real` : the initial approximation p0
- `tol::Float64=1e-5` : the tolerance for error when finding x such that f(x) = 0
- `N::Int64=50` : the max number of iterations to perform
"""
function fixed_point(f::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)
    n = 1

    # while less than max iterations
    while n <= N
        # set current p
        p = f(p0)

        # check difference between p and p0
        if abs(p - p0) < tol
            return p
        end

        # iterate
        n += 1
        # adjust p0 to new p
        p0 = p
    end

    println("Method failed after " * string(N) * "iterations.")
    return nothing
end

"""
    newtons_method(f::Function, g::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)

Find a solution to f(x) = 0 given an intial approximation `p0`.

# Arguments
- `f::Function` : the continuous function given
- `g::Function` : the derivative of the function `f`
- `p0::Real` : the initial approximation p0
- `tol::Float64=1e-5` : the tolerance for error when finding x such that f(x) = 0
- `N::Int64=50` : the max number of iterations to perform
"""
function newtons_method(f::Function, g::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)
    n = 1

    # while less than max iterations
    while n <= N
        # computing p_n
        p = p0 - f(p0) / g(p0)

        # check difference between p and p0
        if abs(p - p0) < tol
            return p
        end

        # iterate
        n += 1
        # adjust p0 to new p
        p0 = p
    end

    println("Method failed after " * string(N) * "iterations.")
    return nothing
end

"""
    secant_method(f::Function, p0::Real, p1::Real; tol::Float64=1e-5, N::Int64=50)

Find a solution to f(x) = 0 given initial approximations `p0` and `p1` such that
x is within the interval [`p0`,`p1`].

# Arguments
- `f::Function` : the continuous function given
- `p0::Real` : the initial approximation p0
- `p1::Real` : the initial approximation p1
- `tol::Float64=1e-5` : the tolerance for error when finding x such that f(x) = 0
- `N::Int64=50` : the max number of iterations to perform
"""
function secant_method(f::Function, p0::Real, p1::Real; tol::Float64=1e-5, N::Int64=50)
    n = 2
    q0 = f(p0)
    q1 = f(p1)

    while n <= N
        p = p1 - (q * (p1 - p0)) / (q1 - q0)

        if abs(p - p1) < tol
            return p
        end

        n += 1

        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
    end

    println("Method failed after " * string(N) * "iterations.")
    return nothing
end

end # end of module
