### A Pluto.jl notebook ###
# v0.19.35

using Markdown
using InteractiveUtils

# ╔═╡ 58fd9e36-892e-11ee-2be1-ef88e2449a09
using DrWatson

# ╔═╡ 0150e304-373d-4dc9-bec9-8b9317fa7f8b
@quickactivate

# ╔═╡ 705257df-ecb3-4a08-8cb5-53b6129ce0e1
begin
	using CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
end
	

# ╔═╡ 3f2f8cd0-5512-4b2f-914c-646a45f94d98
bbfile = ROOTFile(datadir("sims","2nubb_1e8E_EnePhiDist_J23.root"))

# ╔═╡ 446d1348-d0c6-4cfc-93c5-129e410640de
bb = fill_from_root_file(bbfile, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]) 


# ╔═╡ 96492cec-87bf-48e2-8fec-f15f6b8fb6dd
bbSumE       = Process(bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2, "2nubb", :true, SigActivityParams[:bb2Standard], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])


# ╔═╡ 981f4086-5935-45e7-b9a2-71dc74eac9eb
plot(bbSumE.efficiency)

# ╔═╡ faae41c1-1b3c-4b16-942e-0be27bd16430
fig = stephist(bb.reconstructedEnergy2)

# ╔═╡ a2045e81-98c3-4636-8a26-6aa39b89b48e
stephist!(bb.reconstructedEnergy1)

# ╔═╡ 9d5672ea-cfd1-42dd-b6c3-662cac055099
fig

# ╔═╡ 23042b4c-385a-483f-819e-20918c6254dc
default_palettes = Attributes(
    color = Makie.wong_colors(0.1),
    patchcolor = Makie.wong_colors(0.8),
    marker = [:circle, :utriangle, :cross, :rect, :diamond, :dtriangle, :pentagon, :xcross],
    linestyle = [nothing, :dash, :dot, :dashdot, :dashdotdot],
    side = [:left, :right]
)


# ╔═╡ e9051f03-a5c6-4c87-a204-c7d1a96cb126
function dao(basewidth = 2, spinewidth = basewidth * 0.75) 
	dao_attr = Attributes(
        linewidth = basewidth,
		
        Axis = Attributes(
            spinewidth = spinewidth,
            bottomspinevisible = true,
            topspinevisible = true,
            leftspinevisible = true,
            rightspinevisible = true,
            xgridwidth = basewidth,
            ygridwidth = basewidth,
            xgridstyle = :dot,
            ygridstyle = :dot,
			minorgrid= true,
			xminorgridstyle= :dash,
			yminorgridstyle= :dash,
			xautolimitmargin = (0.0f0, 0.0f0),
			yautolimitmargin = (0.0f0, 0.05f0),
        ),
        Legend = Attributes(
            bgcolor = :transparent,
            framevisible = true,
        ),
        Scatter = Attributes(
			linecolor = ["#65ADC2","#111111",  "#233B43", "#E84646", "#C29365", "#362C21", "#316675", "#168E7F", "#109B37"],
            markersize = 15,
            strokewidth = 0,
        ),
		Heatmap = Attributes(
			colormap = :coolwarm,
		),
		palette = default_palettes
		)
    return merge(theme_latexfonts(),dao_attr, Theme(palette = (color=[:red, :blue, :green],)))
end

# ╔═╡ b0e31a6b-b8d1-48ee-8e96-db360cb98952
set_theme!(dao())

# ╔═╡ bca71f87-be13-4e41-899c-224a2da80050
custom_colors = ["red", "blue"]

# ╔═╡ 62540964-0916-4e9c-9408-5ec447f9ba39


# ╔═╡ Cell order:
# ╠═58fd9e36-892e-11ee-2be1-ef88e2449a09
# ╠═0150e304-373d-4dc9-bec9-8b9317fa7f8b
# ╠═705257df-ecb3-4a08-8cb5-53b6129ce0e1
# ╠═3f2f8cd0-5512-4b2f-914c-646a45f94d98
# ╠═446d1348-d0c6-4cfc-93c5-129e410640de
# ╠═96492cec-87bf-48e2-8fec-f15f6b8fb6dd
# ╠═981f4086-5935-45e7-b9a2-71dc74eac9eb
# ╠═faae41c1-1b3c-4b16-942e-0be27bd16430
# ╠═a2045e81-98c3-4636-8a26-6aa39b89b48e
# ╠═9d5672ea-cfd1-42dd-b6c3-662cac055099
# ╠═b0e31a6b-b8d1-48ee-8e96-db360cb98952
# ╠═e9051f03-a5c6-4c87-a204-c7d1a96cb126
# ╠═23042b4c-385a-483f-819e-20918c6254dc
# ╠═bca71f87-be13-4e41-899c-224a2da80050
# ╠═62540964-0916-4e9c-9408-5ec447f9ba39
