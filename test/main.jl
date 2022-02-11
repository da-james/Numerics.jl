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
    h_min = 0.025
    h = 0.05

    u = rkf45(df, u, a, b, h, h_min, 0, tol=0.01)

    display(u)
    plot(u[:,1], u[:,2])
    scatter!(u[:,1], u[:,2])

end


main()
