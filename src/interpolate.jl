"""
    Interpolate

Author: David James, davidabraham@ucla.edu
Date: 20200824
Notes: algorithms and descriptions come from the following
    "Chapter 3: Interpolation and Polynomial Approximation." Numerical Analysis,
        by Richard L. Burden et al., Cengage Learning 2016.

Contains:
- nevilles_method
- newtons_difference
- hermites_method
- cubic_spline
"""
module Interpolate

export nevilles_method

"""
    nevilles_method(x0::Real, x::AbstractVector, q::AbstractVector)

Evaluate the interpolating polynomial P on `n` distinct numbers `x` for the
function `f`.

# Arguments
- `x0::Real` : the value to be approximated
- `x::AbstractVector` : the set of `x` values for the function
- `q::AbstractVector` : the set of f(x) values
"""
function nevilles_method(x0::Real, x::AbstractVector, q::AbstractVector)

    n = size(x)[1]

    for i in n:-1:1
        for j in 1:i
            num = (x0 - x[j]) * q[j+1] - (x0 - x[j+n+1-i]) * q[j]
            dem = x[j+n+1-i] - x[j]
            q[j] = num / dem
        end
    end

    return q[1]
end

"""
    newtons_difference(x::AbstractVector, q::AbstractVector)

Obtain the divided-difference coefficients of the interpolatory polynomial P on
the `n` distinct numbers of `x` for the function `f`.

# Arguments
- `x::AbstractVector` : the set of `x` values for the function
- `q::AbstractVector` : the set of f(x) values
"""
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

"""
    hermites_method(x::AbstractVector, f::AbstractVector, fprime::AbstractVector)

Obtain the coefficients of the Hermite interpolating polynomial H(x) on `n`
distinct numbers of `x` for the function `f`

# Arguments
- `x::AbstractVector` : the set of `x` values for the function
- `f::AbstractVector` : the set of f(x) values
- `fprime::AbstractVector` : the set of f'(x) values
"""
function hermites_method(x::AbstractVector, f::AbstractVector, fprime::AbstractVector)

    n = 2 * size(x)[1]
    z = [x[ceil(Int64, i / 2)] for i in 1:n]
    z1 = [f[ceil(Int64, i / 2)] for i in 1:n]

    q = zeros(n, n)

    # inserting initial values
    q[:,1] = z1

    # inserting first difference values
    z2 = [i%2==1 ? fprime[ceil(Int64, i / 2)] : (q[i+1,1] - q[i,1]) / (z[i+1] - z[i]) for i in 1:n-1]
    pushfirst!(z2, 0)
    q[:,2] = z2

    # creating lower triangular array
    for i in 3:n
        for j in i:n
            q[j,i] = (q[j,i-1] - q[j-1,i-1]) / (z[j] - z[j-i+1])
        end
    end

    return q
end

"""
    cubic_spline(x::AbstractVector, a::AbstractVector)

Construct the cubic spline interpolant `S` for the function `f`, defined within `x`.

# Arguments
- `x::AbstractVector` : the set of `x` values for the function
- `a::AbstractVector` : the set of f(x) values
"""
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

end # end of module
