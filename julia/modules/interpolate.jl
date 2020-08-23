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

    # inserting initial values
    q[:,1] = z1

    # inserting first difference values
    z2 = [i%2==1 ? fprime[ceil(Int64, i / 2)] : (q[i+1,1] - q[i,1]) / (z[i+1] - z[i]) for i in 1:n-1]
    pushfirst!(z2, 0)
    q[:,2] = z2

    # creating lower triangular array
    for i in 3:n
        for j in i:n
            q[j,i] = (q[j,i-1] - q[j-1,i-1]) / (z[j] - z[j-i+1])
        end
    end

    return q
end


end # end of module
