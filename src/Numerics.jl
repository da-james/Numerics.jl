module Numerics
# Author: David James, davidabraham@ucla.edu
# NOTE: each module has references to where the algorithm was
# gathered or learned

include("approximate.jl")
using .Approximate

include("bvp.jl")
using .BVP
export
    finite_difference

include("conversions.jl")
using .Conversions
export
    bar_to_Pa,
    GPa_to_Pa,
    kbar_to_GPa,
    kbar_to_MPa,
    eV_to_J,
    Ry_to_eV,
    A_to_m,
    A_to_Bohr,
    Gyr_to_sec

include("derivative.jl")
using .Derivative
export
    deriv,
    ln_deriv

include("integration.jl")
using .Integration
export
    trapzd,
    simpsn,
    gauleg

include("interpolate.jl")
using .Interpolate
export
    interpolate,
    polint,
    polin2,
    ratint,
    spline,
    splint,
    locate

include("linalg.jl")
using .LinAlg
export
    gauss_elimination,
    crout_factorization,
    gauss_sidel_method,
    norm

include("nde.jl")
using .NDE
export
    jacobian,
    newtons_system

include("ode.jl")
using .ODE
export
    eulers,
    rk4,
    rkf45,
    trapezoid

include("optimization.jl")
using .Optimize
export
    golden_section_search,
    brent_method,
    powell_method

include("pde.jl")
using .PDE

include("polynomials.jl")
using .Polynomials
export
    plgndr,
    psphrc

include("root.jl")
using .Root
export
    bisection_method,
    fixed_point,
    newtons_method,
    secant_method,
    zbrac,
    zbrak

include("util.jl")
using .Util
export
    get_date,
    logspace

end # module
