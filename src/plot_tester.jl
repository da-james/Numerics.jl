using Plots

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i, = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    # linewidth --> range(0, 10, length=n) # makes line thicker
    # seriesalpha --> range(0, 1, length=n) # makes line transparent
    # framestyle --> :none # removes axis
    seriestype --> :scatter
    framestyle --> :origin
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

n = 150
t = range(0, 2π, length=n)
x = sin.(t)
y = cos.(t)

anim = @animate for i ∈ 1:n
    circleplot(x, y, i)
end

gif(anim, "anim_fps15_single.gif", fps=15)
