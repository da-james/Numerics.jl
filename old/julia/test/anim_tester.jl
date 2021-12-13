using Plots

function animation()

    p = plot([sin, cos], zeros(0), leg = false)
    anim = Animation()
    for x = range(0, stop = 10π, length = 100)
        push!(p, x, Float64[sin(x), cos(x)])
        frame(anim)
    end

    gif(anim, "anim_fps30.gif", fps = 30)
end

function animation_2()

    anim = Animation()
    p = plot([sin, cos], 0, π, size=(300, 300))
    scatter!([0], [sin, cos])

    for i in 0:0.1:π
        p[3] = [i], [sin(i)]
        p[4] = [i], [cos(i)]
        frame(anim)
    end

    gif(anim, "anim_fps15.gif", fps=15)

end

animation_2()
