using Plots

using Numerics

function main()

    ode()

end

function ode()

    df(t, y, p) = -2*y[1] + exp(-2*(t-6)^2)
    u = [1]
    a = 0
    b = 10
    h_min = 0.001
    h = 0.1

    u = rkf45(df, u, a, b, h, h_min, 0)
    display(u)

end


main()
