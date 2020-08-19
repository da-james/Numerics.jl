module Interpolate

function nevilles_method(f::Function, x::Real, xarr::Array)
    n = size(xarr)
    q = zeros(n, n)
    q[:,1] = f(xarr)

    for i in 1:n
        for j in 1:i
            num = (x - xarr[i-j+1]) * q[i+1,j+1] - (x - xarr[i+1]) * q[i,j]
            dem = x - xarr[i-j+1]
            q[i+1,j+1] = num / dem
        end
    end

    return q
end


end # end of module
