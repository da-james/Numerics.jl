module Derivative

export deriv

function deriv(f, x0::Real, p, h::Real=0.001)

    k1 = f(x0 - 2h, p...)
    k2 = f(x0 - h, p...)
    k3 = f(x0 + h, p...)
    k4 = f(x0 + 2h, p...)

    return 1/(12*h) * (k1 - 8*k2 + 8*k3 - k4)
end

end # module
