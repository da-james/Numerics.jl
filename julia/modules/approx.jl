"""
    Approximate

Author: David James, davidabraham@ucla.edu
Date: 20200829
Notes: algorithms and descriptions come from the following
    "Chapter 8: Approximation Theory." Numerical Analysis,
        by Richard L. Burden et al., Cengage Learning 2016.

Contains:
- fft
"""
module Approximate

"""
    fft(x::AbstractVector)

Compute the coefficients in the summation
    1/m*Σ{k=0, 2*m-1}{c_k * exp((k*x)im)}
"""
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

end # end of module
