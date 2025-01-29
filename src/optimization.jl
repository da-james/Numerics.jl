module Optimize

using ..Numerics: norm

export golden_section_serach, brent_method, powell_method

# Golden section constants
const PHI = (1 + √5) / 2
const RESPHI = 2 - PHI  # 1/phi

# test functions for brent when I properly incorporate testing (lul)
# easier
# f(x) = (x - 2)^2 + 3; [-5, 5]; x_min = 2
# more complex with a local minimum
# g(x) = x^4 - 3x^3 + 2; [-1, 3]; x_min = 2.25

# function quadratic(x)
#     return (x[1] - 3)^2 + (x[2] + 2)^2
# end
# x0 = [0.0, 0.0] --> [3, -2]

# function rosenbrock(x)
#     return 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
# end
# x0 = [-1.2, 1.0] --> [1.0, 1.0]

# function sum_of_squares(x)
#     return x[1]^2 + x[2]^2
# end
# x0 = [1.0, 1.0] --> [0.0, 0.0]

# function hartmann(x)
#     return -((1.1 - 3*x[1] + x[1]^2)^2 + (1.1 - 3*x[2] + x[2]^2)^2)
# end
# x0 = [0.5, 0.5] --> [1.5, 1.5]



"""
    golden_section_search(f, a, b; tol=1e-5)

Finds the minimum of `f` in the interval `[a, b]` using the golden section search.

"""
function golden_section_search(f, a, b; tol=1e-5)
    x1 = a + RESPHI * (b - a)
    x2 = b - RESPHI * (b - a)
    f1, f2 = f(x1), f(x2)

    while abs(b - a) > tol
        if f1 < f2
            b, x2, f2 = x2, x1, f1
            x1 = a + RESPHI * (b - a)
            f1 = f(x1)
        else
            a, x1, f1 = x1, x2, f2
            x2 = b - RESPHI * (b - a)
            f2 = f(x2)
        end
    end

    return (a + b) / 2
end

"""
    brent_method(f, a, b; tol=1e-5)

Finds the minimum of `f` in the interval `[a, b]` using Brent's method.
"""
function brent_method(f, a, b; tol=1e-5)
    golden_ratio = PHI - 1
    x = w = v = (a + b) / 2
    fx = fw = fv = f(x)
    d = e = b - a

    while abs(b - a) > tol
        g = e
        e = d
        u = nothing

        # Attempt parabolic fit
        if x != w && x != v && w != v
            u = x - ((x - w)^2 * (fx - fv) - (x - v)^2 * (fx - fw)) /
                  (2 * ((x - w) * (fx - fv) - (x - v) * (fx - fw)))
            if a + tol <= u <= b - tol && abs(u - x) < g / 2
                d = abs(u - x)
            else
                u = nothing
            end
        end

        if u === nothing
            if x < (a + b) / 2
                u = x + golden_ratio * (b - x)
                d = b - x
            else
                u = x - golden_ratio * (x - a)
                d = x - a
            end
        end

        fu = f(u)
        if fu <= fx
            if u < x
                b = x
            else
                a = x
            end
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else
            if u < x
                a = u
            else
                b = u
            end
            if fu <= fw || w == x
                v, fv = w, fw
                w, fw = u, fu
            elseif fu <= fv || v == x || v == w
                v, fv = u, fu
            end
        end
    end

    return x
end

"""
    powell_method(f, x0; directions=nothing, tol=1e-6, max_iters=100)

Performs Powell's method for multidimensional minimization of `f` starting at `x0`.
If `directions` is not provided, the identity matrix is used as the initial directions.
"""
function powell_method(f, x0; directions=nothing, tol=1e-6, max_iters=100)
    n = length(x0)
    # Use provided directions or default to identity matrix
    if directions === nothing
        directions = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    else
        # Ensure the provided directions matrix has the correct dimensions
        if size(directions) != (n, n)
            error("Provided directions must be an $n x $n matrix.")
        end
    end

    x = copy(x0)
    fx = f(x)

    for iter in 1:max_iters
        x_start = copy(x)
        fx_start = fx

        # Iterate through all directions
        for i in 1:n
            direction = directions[:, i]

            # Dynamic line search range
            f_line = α -> f(x + α * direction)
            α_min = brent_method(f_line, -10.0, 10.0; tol=tol)
            x .= x + α_min * direction
            fx = f(x)
        end

        # Check for convergence
        if norm(x - x_start) / max(norm(x_start), tol) < tol &&
           abs(fx - fx_start) / max(abs(fx_start), tol) < tol
            println("Converged at iteration $iter with x = $x, f(x) = $fx")
            return x, fx
        end

        # Update directions using new composite direction
        new_direction = x - x_start
        if norm(new_direction) > tol
            new_direction ./= norm(new_direction)  # Normalize
        else
            new_direction .= directions[:, end]  # Use last direction
        end

        # Shift and rescale directions
        directions[:, 1:end-1] .= directions[:, 2:end]
        directions[:, end] = new_direction
        for i in 1:n
            directions[:, i] ./= norm(directions[:, i])  # Normalize directions
        end
    end

    println("Maximum iterations reached. Current minimum at x = $x, f(x) = $fx")
    return x, fx
end


end  # module
