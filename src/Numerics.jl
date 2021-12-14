module Numerics

include("approx.jl")
include("bvp.jl")
include("interpolate.jl")
include("nde.jl")
include("ode.jl")
include("pde.jl")
include("roots.jl")

using Approximate, BVP, Interpolate, LinAlg, NDE, ODE, PDE, Roots

export
    finite_difference,
    nevilles_method,
    gauss_elimination,
    crout_factorization,
    gauss_sidel_method,
    norm,
    jacobian,
    newtons_system,
    rk4_method,
    rk4_system,
    trapezoid_method,
    biscetion_method,
    fixed_point,
    newtons_method,
    secant_method

end # module
