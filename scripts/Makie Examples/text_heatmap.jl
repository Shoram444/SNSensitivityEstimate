using CairoMakie, Random

Random.seed!(123)
m = 11
n = 11
data = rand(m, n)

ticks = string.(0:100:1000)

fig = Figure(
    size = (600, 600),
    fontsize = 18,
    pt_per_unit = 1
)

ax = Axis(
    fig[1,1],
    xticks = (1:m, ticks),
    yticks = (1:n, ticks),
    xlabel = L"E_{min} ~[keV]",
    ylabel = L"E_{max} ~[keV]",
    xlabelsize = 26,
    ylabelsize = 26,
    aspect = 1,
    title = "title"
)
hmap = heatmap!(ax, data, colormap =:plasma)
for i in 1:11, j in 1:11
    txtcolor = data[i, j] < 0.15 ? :white : :black
    text!(
        ax, 
        "$(round(data[i,j], digits = 2))", 
        position = (i, j),
        color = txtcolor, 
        align = (:center, :center)
    )
end
Colorbar(fig[1, 2], hmap; label = "value", width = 15, ticksize = 15)
colsize!(fig.layout, 1, Aspect(1, 1))
resize_to_layout!(fig)
fig
save("fig.svg", fig, pt_per_unit = 1)

with_theme(theme_latexfonts()) do 
    fig
    Label(fig[0,1], "title label", tellwidth = false)
    fig
end

