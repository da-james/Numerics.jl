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

function newtons_difference(x::AbstractArray, q::AbstractArray)

    n = size(x)[1]
    f = zeros(n, n)
    f[:,1] = q

    for i in 2:n
        for j in i:n
            f[j,i] = (f[j,i-1] - f[j-1,i-1]) / (x[j] - x[j-i+1])
        end
    end

    return f
end

function hermites_method(x::AbstractArray, f::AbstractArray, fprime::AbstractArray)

    n = 2 * size(x)[1]
    z = [x[ceil(Int64, i / 2)] for i in 1:n]
    z1 = [f[ceil(Int64, i / 2)] for i in 1:n]

    q = zeros(n, n)
    q[:,1] = z1

    for r in 1:n-2
        for c in 2:n
            if c==2 && r%2==1
                q[r,c] = fprime[ceil(Int64, r / 2)]
            elseif c==2 && r%2==0
                q[r,c] = (q[r+1,c-1] - q[r,c-1]) / (z[r+1] - z[r])
            else
                q[r,c] = (q[r+1,c-1] - q[r,c-1]) / (z[r+2] - z[r])
            end
        end
    end

    for row in eachrow(q)
        println(row)
    end

end


end # end of module
