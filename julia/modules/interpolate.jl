module Interpolate

function nevilles_method(f::Function, x::Real, xarr::Array)
    n = size(xarr)
    q = zeros(n, n)

    for i in 1:n
        for j in 1:i
            num = (x - xarr[i-j]) * q[i,j-1] - (x - xarr[i]) * q[i-1,j-1]
            dem = x - xarr[i-j]
            q[i,j] = num / dem
        end
    end

    return q
end


end # end of module
