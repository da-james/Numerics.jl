
module Approximate

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
