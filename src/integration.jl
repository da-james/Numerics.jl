module Integration

export trapzd, simpsn, gauleg

"""
    trapzd(f::Function, x1::Real, x2::Real, n::Int)

Integration method using the trapezoidal rule.
"""
function trapzd(f::Function, x1::Real, x2::Real, n::Int)

    h = (x2 - x1) / n

    sum = 0
    for i in 1:n
        sum += h/2 * (f(x1) + f(x2))
    end

    return sum
end

"""
    simpsn(f::Function, x1::Real, x2::Real, n::Int)

Integration method using Simpson's rule.
"""
function simpsn(f::Function, x1::Real, x2::Real, n::Int)

    h = (x2 - x1) / n
    xh = x1 + h

    sum = 0
    for i in 1:n
        sum += h/3 * (f(x1) + 4*f(xh) + f(x2))
    end

    return sum
end


"""
    gauleg(x1::Real, x2::Real, n::Int, reltol::Real=3e-14)

Given the lower and upper limits of integration `x1` and `x2`, and
given `n`, this routine returns an array of abscissas and weights of
the Gauss-Legendre n-point quadrature formula

# Arguments
- `x1::Real`           : lower limit of integration
- `x2::Real`           : upper limit of integration
- `n::Int`             : length of array to return
- `reltol::Real=3e-14` : tolerance to solve out weights
"""
function gauleg(x1::Real, x2::Real, n::Int, reltol::Real=3e-14)

    m = floor(Int64, (n + 1) / 2)
    xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)

    x = zeros(n)
    w = zeros(n)

    for i in 1:m

        z = cos(π * (i - 0.25) / (n + 0.5))
        pp = 0

        while true

            p1 = 1
            p2 = 0

            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2*j - 1)*z*p2 - (j - 1)*p3) / j
            end

            pp = n*(z*p1 - p2) / (z^2 - 1)
            z1 = z
            z = z1 - p1 / pp

            (abs(z - z1) > reltol) || break
        end

        x[i] = xm - xl*z
        x[n+1-i] = xm + xl*z

        w[i] = 2*xl / ((1 - z^2)*pp^2)
        w[n+1-i] = w[i]

    end

    return [x w]
end

end # module
