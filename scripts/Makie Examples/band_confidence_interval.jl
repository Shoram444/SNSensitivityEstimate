using CairoMakie

x = y = -10:0.11:10
y1d = sin.(x) ./ x

lower = y1d .- 0.1
upper = y1d .+ 0.1

fig = Figure(size= (800, 600))
ax = Axis(
    fig[1,1],
    xlabel = "x",
    ylabel = "y",
    title = "1/sin(x) ± 0.1",
)
lines!(ax, x, y1d, color = :blue)
band!(ax, x , lower, upper, color =(:grey, 0.2))
xlims!(extrema(x)...)
fig