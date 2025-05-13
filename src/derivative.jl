module Derivative

export deriv

function deriv(f, x0::Real, args...; h::Real=1e-3)

    k1 = f(x0 - 2h, args...)
    k2 = f(x0 - h, args...)
    k3 = f(x0 + h, args...)
    k4 = f(x0 + 2h, args...)

    return (k1 - 8*k2 + 8*k3 - k4) / (12*h)
end

end # module
