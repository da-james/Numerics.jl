module ODE

function eulers_method(f::Function, α::Real, a::Real, b::Real, n::Int64)

    n1 = n+1
    u = zeros(n1,2)

    h = (b - a) / n
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        u[i,2] = u[i-1,2] + h * f(u[i-1,1], u[i-1,2])
        u[i,1] = a + (i - 1) * h
    end

    return u
end



end # end of module
