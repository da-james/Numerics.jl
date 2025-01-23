using Plots
include("Numerics.jl")

function quadratic(x)
    return (x[1] - 3)^2 + (x[2] + 2)^2
end

function rosenbrock(x)
    return 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
end

function main()

    x0 = [-1.2, 1.0]
    x1, f = Numerics.Optimize.powell_method(quadratic, x0)
    println(x1)
    println(f)

    xs = collect(-2.0:0.05:4.0)
    ys = collect(-4:0.1:3.0)
    X = [x for x = xs for _ = ys]
    Y = [y for _ = xs for y = ys]
    V = [X Y]
    Z = zeros(length(V[:,1]))
    for i in 1:length(V[:,1])
        Z[i] = quadratic([V[i,1], V[i,2]])
    end


    plt = surface(X, Y, Z, xlabel = "x",
                  ylabel = "y",
                  zlabel = "z", camera = (10, round(atand(1 / √2); digits=3)), colorbar=:none,  dpi=600)
    plt = scatter!([x0[1]], [x0[2]], [quadratic(x0)], label="x0")
    plt = scatter!([x1[1]], [x1[2]], [f], label="minimum")
    display(plt)

end
