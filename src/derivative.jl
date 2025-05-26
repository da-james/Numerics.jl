module Derivative

export deriv, ln_deriv

function deriv(f, x0::Real, args...; h::Real=1e-3, rel::Real=1e-3)

    h₀ = max(h, rel*abs(x0))
    k1 = f(x0 - 2h₀, args...)
    k2 = f(x0 - h₀, args...)
    k3 = f(x0 + h₀, args...)
    k4 = f(x0 + 2h₀, args...)

    return (k1 - 8*k2 + 8*k3 - k4) / (12*h₀)
end

function ln_deriv(f, x0, args...; h::Real=1e-3, rel::Real=1e-3)
    return deriv(x -> log(f(x, args...)), x0; h=h, rel=rel)
end

end # module
