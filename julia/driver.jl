include("modules/roots.jl")

using .Root


function main()

    f(x) = x - (x^3 + 4*x^2 - 10) / (3*x^2 + 8*x)
    # zero = Root.bisection_method(f, -2, -3)
    sol = Root.fixed_point(f, 1.5)
    println(sol)

end

main()
