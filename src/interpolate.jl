"""
    Interpolate

Author: David James, davidabraham@ucla.edu
Date: 20211221
Notes: algorithms and descriptions come from the following
    "Chapter 3: Interpolation and Extrapolation" Numerical Recipes in
    FORTRAN 77 by Press W. H. et al.

Contains:
- polint
- ratint
- spline
- splint
"""
module Interpolate

export nevilles_method

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

        drift = abs(x - xa[i])
        if drift < diff
            ns = i
            diff = drift
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
            c[i] = hp * den
            d[i] = ho * den

        end

        # after each column in the tableau is completed, we decide
        # which correction, c or d, we want to add to our accumulating
        # value y, i.e. which path to take through the tableau
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


end # end of module
