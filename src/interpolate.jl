"""
    Interpolate

Author: David James, davidabraham@ucla.edu
Date: 20220627
Notes: algorithms and descriptions come from the following
    "Chapter 3: Interpolation and Extrapolation" Numerical Recipes in
    FORTRAN 77 by Press W. H. et al.

Contains:
- interpolate
- polint
- ratint
- spline
- splint
"""
module Interpolate

export interpolate, polint, ratint, spline, splint, locate

"""
    interpolate(xa::AbstractVector, ya::AbstractVector, x::Real; k::Int=2)

Wrapper for `polint` function. This locates the index of the wanted `x` in `xa`
and interpolates the `ya` based on the desired degree, `k`, of interpolation.
For more information refer to `polint` and/or `locate` for the operations
performed.
"""
function interpolate(xa::AbstractVector, ya::AbstractVector, x::Real; k::Int=2)

    i = locate(xa, x)

    if(i == 1)
        y, dy = polint(xa[i:i+k], ya[i:i+k], x)
    elseif(i == length(xa) || i == length(xa) - k + 1)
        y, dy = polint(xa[i-k:i], ya[i-k:i], x)
    else
        l = k - 1
        y, dy = polint(xa[i-l:i+l], ya[i-l:i+l], x)
    end

    return (y, dy)
end

"""
    polint(xa::AbstractVector, ya::AbstractVector, x::Real)

Given arrays `xa` and `ya` each of length `n`, and given a value `x`,
this routine returns a value `y`, and an error estimate `dy`. If `P(x)`,
is a polynomial of degree `N - 1` such that `P(xaᵢ) = yaᵢ, i = 1,...,n`,
then the returned value `y = P(x)`.
"""
function polint(xa::AbstractVector, ya::AbstractVector, x::Real)

    n = length(xa)
    c = zeros(n)
    d = zeros(n)

    ns = 1
    diff = abs(x - xa[1])
    dy = 0

    # here we find the index ns of the closest table entry
    for i in 1:n

        dift = abs(x - xa[i])
        if dift < diff
            ns = i
            diff = dift
        end

        # and initialize the tableau of c's and d's
        c[i] = ya[i]
        d[i] = ya[i]

    end

    # this is the initial approximation of y
    y = ya[ns]
    ns -= 1

    for m in 1:n-1
        # for each column of the tableau
        # we loop over the current c's and d's and update them
        for i in 1:n-m

            ho = xa[i] - x
            hp = xa[i+m] - x
            den = ho - hp

            w = c[i+1] - d[i]

            # if den == 0 --> exit out

            den = w / den

            # c's and d's get updated
            c[i] = ho * den
            d[i] = hp * den

        end

        # after each column in the tableau is completed, we decide
        # which correction, c or d, we want to add to our accumulating
        # value y, i.e. which path to take through the tableau
        if 2*ns < n - m
            dy = c[ns+1]
        else
            dy = d[ns]
            ns -= 1
        end

        y += dy

    end

    return (y, dy)
end

"""
    ratint(xa::AbstrractVector, ya::AbstractVector, x::Real)

Given arrays of `xa` and `ya` each of length `n`, and given a value `x`,
this routine returns a value of `y` and an accuracy estimate of `dy`.
The value returned is that of the diagonal rational function,
evaluated at `x`, which passes through the `n` points `(xaᵢ, yaᵢ),i=1,...,n`.
"""
function ratint(xa::AbstractVector, ya::AbstractVector, x::Real)

    n = length(xa)
    c = zeros(n)
    d = zeros(n)
    ϵ = 1e-25

    ns = 1
    hh = abs(x - xa[1])
    dy = 0

    for i in 1:n

        h = abs(x - xa[i])

        if h == 0
            y = ya[i]
            dy = 0
            return (y, dy)
        elseif h < hh
            ns = i
            hh = h
        end

        c[i] = ya[i]
        d[i] = ya[i] + ϵ

    end

    y = ya[ns]
    ns = ns - 1

    for m in 1:n-1
        for i in 1:n-m

            w = c[i+1] - d[i]
            # h will never be zero, since this was tested in the
            # initializing loop
            h = xa[i+m] - x
            t = (xa[i] - x) * d[i] / h

            dd = t - c[i+1]
            # if dd == 0 --> exit out

            dd = w / dd

            c[i] = t * dd
            d[i] = c[i+1] * dd

        end

        if 2*ns < n - m
            dy = c[ns+1]
        else
            dy = d[ns]
            ns = ns - 1
        end

        y += dy

    end

    return (y, dy)
end

"""
    spline(x::AbstractVector, y::AbstractVector, yp1::Real, yp2::Real)

Given arrays `x[1:n]` and `y[1:n]` containing a tabulated function,
i.e. `yᵢ=f(xᵢ)`, with `x₁<x₂<...<xₙ`, and given `yp1` and `ypn` for
the first derivative of the interpolating function at points 1 and `n`,
respectively, this routine returns an array `y2[1:n]` of length `n`
which contains the second derivatives of the interpolating function
at the tabulated points `xᵢ`. If `yp1` or `ypn` are equal to 1e30 or
larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that
boundary.
"""
function spline(x::AbstractVector, y::AbstractVector, yp1::Real, yp2::Real)

    n = length(x)
    u = zeros(n)
    y2 = zeros(n)

    if yp1 > 99e30
        y2[1] = 0
        u[1] = 0
    else
        y2[1] = -0.5
        u[1] = (3 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1)
    end

    for i in 2:n-1

        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p = sig * y2[i-1] + 2

        y2[i] = (sig - 1) / p

        a = (y[i+1] - y[i]) / (x[i+1] - x[i])
        b = (y[i] - y[i-1]) / (x[i] - x[i-1])

        u[i] = 6*(a - b) / (x[i+1] - x[i-1]) - sig*u[i-1] / p

    end

    if ypn > 99e30
        qn = 0
        un = 0
    else
        qn = 0.5

        a = 3 / (x[n] - x[n-1])
        b = (ypn - (y[n] - y[n-1]) / (x[n] - x[n-1]))
        un = a * b
    end

    y2[n] = (un - qn*u[n-1]) / (qn*y2[n-1] + 1)

    for i in n-1:-1:1

        y2[i] = y2[i] * y2[i+1] + u[i]

    end

    return y2

end

"""
    splint(xa::AbstractVector, ya::AbstractVector, y2a::AbstractVector, x::Real)

Given the arrays `xa[1:n]` and `ya[1:n]` of length `n`, which tabulate
a function (with the xaᵢ's in order), and given the array `y2a[1:n]`,
which is the output from the `spline` function, and given a value `x`,
this routine returns a cubic-spline interpolated value `y`
"""
function splint(xa::AbstractVector, ya::AbstractVector, y2a::AbstractVector, x::Real)

    n = length(xa)

    klo = 1
    khi = n

    while khi - klo > 1

        k = (khi + klo) / 2
        if xa[k] > x
            khi = k
        else
            klo = k
        end

    end

    h = xa[khi] - xa[klo]
    # if h == 0 --> exit out

    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    c = ((a^3 - a)*y2a[klo] + (b^3 - b)*y2a[khi]) * h^2 / 6
    y = a*ya[klo] + b*ya[khi] + c

    return y

end

"""
    locate(data::AbstractVector, mark::Real)

Finds the index location of the given `mark` in the vector `data`.
If the `mark` is in-between values in `data` then the location given
is the left bound s.t. `a < mark < b` where `a` is the left bound.
If the value is not in the then the location given is either the
minimum or maximum of `data` depending on whether `mark` is larger
relatively to the elements inside `data`. This function should be
used in conjunction with one of the interpolation calls to have an `x`
to insert into `polint` or other scheme in the library.
"""
function locate(data::AbstractVector, mark::Real)

    diff = abs.(data .- mark)
    val, loc = findmin(diff)

    return loc

end

end # end of module
