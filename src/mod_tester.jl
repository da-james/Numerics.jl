include("modules/roots.jl")
include("modules/interpolate.jl")
include("modules/linalg.jl")
include("modules/ode.jl")
include("modules/approx.jl")
include("modules/nonlinear.jl")
include("modules/bvp.jl")
include("modules/pde.jl")

using .Root
using .Interpolate
using .LinAlg
using .ODE
using .Approximate
using .Nonlinear
using .BVP
using .PDE


function math151a()

    # f(x) = x - (x^3 + 4*x^2 - 10) / (3*x^2 + 8*x)
    # zero = Root.bisection_method(f, -2, -3)
    # sol = Root.fixed_point(f, 1.5)
    # println(sol)

    # x = 1940:10:1990
    # q = [132165, 151326, 179323, 203302, 226542, 249633]
    # println(Interpolate.nevilles_method(2000, 5, x, q))

    # x = [8.1, 8.3, 8.6, 8.7]
    # y = [16.9446, 17.56492, 18.50515, 18.82091]
    # println(Interpolate.newtons_difference(x, y))

    # x = [1.3, 1.6, 1.9]
    # f = [0.6200860, 0.4554022, 0.2818186]
    # fprime = [-0.5220232, -0.5698959, -0.5811571]
    # println(Interpolate.hermites_method(x, f, fprime))

    # x = [0.9, 1.3, 1.9, 2.1, 2.6, 3.0, 3.9, 4.4, 4.7, 5.0, 6.0,
    #      7.0, 8.0, 9.2, 10.5, 11.3, 11.6, 12.0, 12.6, 13.0, 13.3]
    # a = [1.3, 1.5, 1.85, 2.1, 2.6, 2.7, 2.4, 2.15, 2.05, 2.1, 2.25,
    #      2.3, 2.25, 1.95, 1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25]
    # println(Interpolate.cubic_spline(x, a))

    # a = [1 -1 2 -1 -8;
    #      2 -2 3 -3 -20;
    #      1 1 1 0 -2;
    #      1 -1 4 3 4]
    # println(LinAlg.gauss_elimination(a))

    # a = [2 -1 0 0 1;
    #      -1 2 -1 0 0;
    #      0 -1 2 -1 0;
    #      0 0 -1 2 1]
    # println(LinAlg.crout_factorization(a))

    a = [10 -1 2 0;
         -1 11 -1 3;
         2 -1 10 -1;
         0 3 -1 8]
    b = [6, 25, -11, 15]
    x0 = [0, 0, 0, 0]
    println(LinAlg.gauss_sidel_method(a,b,x0))
end

function math151b()
    # df(t, y) = y - t^2 + 1
    # a = 0
    # b = 2
    # y0 = 0.5
    # n = 10
    # println(ODE.eulers_method(df, y0, a, b, n))
    # arr = ODE.rk4_method(df, y0, a, b, n)
    # display(arr)

    # function dI(t, i)
    #     i1 = i[1]
    #     i2 = i[2]
    #     di1 = -4 * i1 + 3 * i2 + 6
    #     di2 = -2.4 * i1 + 1.6 * i2 + 3.6
    #     return [di1, di2]
    # end
    # a = 0
    # b = 0.5
    # i0 = [0, 0]
    # n = 5
    # display(ODE.rk4_system(dI, i0, a, b, n))

    # f(t, y) = 5 * exp(5 * t) * (y - t)^2 + 1
    # df(t, y) = 10 * exp(5 * t) * (y - t)
    # a = 0
    # b = 1
    # y0 = -1
    # n = 4
    # display(ODE.trapezoid_method(f, df, y0, a, b, n, tol=1e-6))
    # display(ODE.rk4_method(f, y0, a, b, n))

    # y = [0.26440, 0.84081, 1.36150, 1.61282, 1.36672, 0.71697, 0.07909, -0.14576]
    # display(Approximate.fft(y))

    # function df(x)
    #     x1 = x[1]
    #     x2 = x[2]
    #     x3 = x[3]
    #     df1 = 3 * x1 - cos(x2 * x3) - 1 / 2
    #     df2 = x1^2 - (81 * (x2 + 0.1)^2) + sin(x3) + 1.06
    #     df3 = exp(-x1 * x2) + 20 * x3 + ((10 * pi - 3) / 3)
    #     return [df1, df2, df3]
    # end
    # x0 = [0.1, 0.1, -0.1]
    # display(Nonlinear.newtons_system(df, x0))

    # p(x) = -2 / x
    # q(x) = 2 / x^2
    # r(x) = sin(log(x)) / x^2
    # a = 1
    # b = 2
    # α = 1
    # β = 2
    # display(BVP.finite_difference(p, q, r, a, b, α, β, 9))

    x = (0, 2)
    y = (0, 1)
    f(x,y) = x
    g(x,y) = exp(x)
    n = 6
    m = 5
    display(PDE.poisson_finite_difference(f, g, x, y, m, n))

end

math151b()
