module LinAlg

function gauss_elimination(a::AbstractArray)

    n = size(a)[1]

    for i in 1:n-1
        p = min(i, n)
        c = i
        while a[p,i] == 0
            c += 1
            p = min(c, n)
            if c > n
                println("no unique solution exists")
                return
            end
        end

        if p != i
            r = a[p,:]
            a[p,:] = a[i,:]
            a[i,:] = r
        else
            for j in i+1:n
                m = a[j,i] / a[i,i]
                r = a[j,:] - m .* a[i,:]
                a[j,:] = r
            end
        end
    end

    if a[n,n] == 0
        println("no unique solution exists")
        return
    end

    x = zeros(n)
    x[n] = a[n,n+1] / a[n,n]

    for i in n-1:-1:1
        xj = 0
        for j = i+1:n
            xj += a[i,j] * x[j]
        end

        x[i] = (a[i,n+1] - xj) / a[i,i]
    end

    return x
end

function crout_factorization(a::AbstractArray)
    n = size(a)[1]

    l = zeros(n, n)
    u = zeros(n, n)
    z = zeros(n)
    x = zeros(n)

    l[1,1] = a[1,1]
    u[1,2] = a[1,2] / l[1,1]
    z[1] = a[1,n+1] / l[1,1]

    for i in 2:n-1
        l[i,i-1] = a[i,i-1]
        l[i,i] = a[i,i] - l[i,i-1] * u[i-1,i]
        u[i,i+1] = a[i,i+1] / l[i,i]
        z[i] = (a[i,n+1] - l[i,i-1] * z[i-1]) / l[i,i]
    end

    l[n,n-1] = a[n,n-1]
    l[n,n] = a[n,n] - l[n,n-1] * u[n-1,n]
    z[n] = (a[n,n+1] - l[n,n-1] * z[n-1]) / l[n,n]

    x[n] = z[n]

    for i in n-1:-1:1
        x[i] = z[i] - u[i,i+1] * x[i+1]
    end

    return x
end

function gauss_sidel_method(a::AbstractArray, b::AbstractArray, x0::AbstractArray; tol::Float64=1e-5, N::Int64=50)

    n = size(a)[1]
    x = zeros(n)

    k = 1
    while k <= N
        for i in 1:n
            first = 0
            for j in 1:i-1
                first += a[i,j] * x[j]
            end
            second = 0
            for j in i+1:n
                second += a[i,j] * x0[j]
            end

            x[i] = 1 / a[i,i] * (-1 * first -1 * second + b[i])
        end

        if norm(x, x0) < tol
            return x
        end

        k += 1
        x0 = x

    end

    println("Maximum number of iterations exceeded")
    return
end

function norm(x::AbstractArray, y::AbstractArray)

    n = size(x)[1]

    total = 0
    for i in 1:n
        total += (x[i] - y[i])^2
    end

    return sqrt(total)
end


end # end of module
