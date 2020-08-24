include("modules/roots.jl")
include("modules/interpolate.jl")
include("modules/linalg.jl")

using .Root
using .Interpolate
using .LinAlg


function main()

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

main()
