module ODE

function eulers_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

    n1 = N + 1
    u = zeros(n1,2)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        u[i,2] = u[i-1,2] + h * f(u[i-1,1], u[i-1,2])
        u[i,1] = a + (i - 1) * h
    end

    return u
end

function rk4_method(f::Function, α::Real, a::Real, b::Real, N::Int64)

    n1 = N + 1
    u = zeros(n1,2)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        t = u[i-1,1]
        w = u[i-1,2]

        k1 = h * f(t, w)
        k2 = h * f(t + h / 2, w + k1 / 2)
        k3 = h * f(t + h / 2, w + k2 / 2)
        k4 = h * f(t + h, w + k3)

        u[i,2] = w + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        u[i,1] = a + (i - 1) * h
    end

    return u
end

function rk4_system(f::Function, α::AbstractArray, a::Real, b::Real, N::Int64)

    n1 = N + 1
    m = size(α)[1]

    u = zeros(n1, m + 1)

    h = (b - a) / N
    u[1,1] = a
    u[1,2:end] = α

    for i in 2:n1
        t = u[i-1,1]

        v = u[i-1,2:end]

        k1 = h * f(t, v)
        k2 = h * f(t + h / 2, v .+ k1 / 2)
        k3 = h * f(t + h / 2, v .+ k2 / 2)
        k4 = h * f(t + h, v .+ k3)

        u[i,2:end] = v .+ (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) ./ 6
        u[i,1] = a + (i - 1) * h
    end

    return u
end

function trapezoid_method(f::Function, df::Function, α::Real, a::Real, b::Real, N::Int64; tol::Float64=1e-5, M::Int64=10)

    n1 = N + 1
    u = zeros(n1,2)

    h = (b - a) / N
    u[1,1] = a
    u[1,2] = α

    for i in 2:n1
        ti = u[i-1,1]
        wi = u[i-1,2]

        k1 = wi + h / 2 * f(ti, wi)
        w0 = k1
        j = 1
        FLAG = 0

        while FLAG == 0
            num = w0 - (h / 2 * f(ti + h, w0)) - k1
            dem = 1 - h / 2 * df(ti + h, w0)

            w = w0 - num / dem

            if abs(w - w0) < tol
                FLAG = 1
            else
                j += 1
                w0 = w
                if j > M
                    println("The maximum number of iterations exceeded")
                    return nothing
                end
            end
            u[i, 2] = w
        end

        u[i, 1] = a + (i - 1) * h
    end

    return u
end

end # end of module
