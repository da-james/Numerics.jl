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

end # end of module
