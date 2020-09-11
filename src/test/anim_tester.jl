using Plots

p = plot([sin, cos], zeros(0), leg = false)
anim = Animation()
for x = range(0, stop = 10π, length = 100)
    push!(p, x, Float64[sin(x), cos(x)])
    frame(anim)
end

gif(anim, "anim_fps30.gif", fps = 30)
