module Interpolate

function nevilles_method(x0::Real, n::Int64, x::AbstractArray, q::Array)

    for i in n:-1:1
        for j in 1:i
            num = (x0 - x[j]) * q[j+1] - (x0 - x[j+n+1-i]) * q[j]
            dem = x[j+n+1-i] - x[j]
            q[j] = num / dem
        end
    end

    return q[1]
end


end # end of module
