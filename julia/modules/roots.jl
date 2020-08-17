module Root

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
            println("iterations: ", n)
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
    return -1
end

function fixed_point(f::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)
    n = 1

    while n <= N
        p = f(p0)

        if abs(p - p0) < tol
            return p
        end

        n += 1
        p0 = p
    end

    println("Method failed after " * string(N) * "iterations.")
    return -1
end

function newtons_method(f::Function, g::Function, p0::Real; tol::Float64=1e-5, N::Int64=50)
    n = 1

    while n <= N
        p = p0 - f(p0) / g(p0)

        if abs(p - p0) < tol
            return p
        end

        n += 1
        p0 = p
    end

    println("Method failed after " * string(N) * "iterations.")
    return -1
end

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
    return -1
end


end # end of module
