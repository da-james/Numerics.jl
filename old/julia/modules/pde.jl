module PDE

function poisson_finite_difference(f::Function, g::Function, x::Tuple, y::Tuple, m::Int64, n::Int64; tol::Float64=1e-5, N::Int64=50)

    h = (x[2] - x[1]) / n
    k = (y[2] - y[1]) / m

    x = (x[1] + h:h:x[2] - h)
    y = (y[1] + k:k:y[2] - k)

    for i in n
        for j in m
            l = i + (m - 1 - j) * (n - 1)

        end
    end


end

end # end of module
