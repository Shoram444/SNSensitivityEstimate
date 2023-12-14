using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

Bi214file = ROOTFile("data/sims/Bi214_EnePhiDistPos_J23.root")
Tl208file = ROOTFile("data/sims/Tl208_EnePhiDistPos_J23.root")
K40file = ROOTFile("data/sims/K40_EnePhiDistPos_J23.root")
Pa234mfile = ROOTFile("data/sims/Pa234m_EnePhiDistPos_J23.root")
bbfile = ROOTFile("data/sims/2nubb_EnePhiDistPos_J23.root")   # Standard 2nubb


Bi214 = fill_from_root_file(Bi214file, "tree", keys(Bi214file["tree"]))     |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
Tl208 = fill_from_root_file(Tl208file, "tree", keys(Tl208file["tree"]))     |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
K40 = fill_from_root_file(K40file, "tree", keys(K40file["tree"]))           |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
Pa234m = fill_from_root_file(Pa234mfile, "tree", keys(Pa234mfile["tree"]))  |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
bb = fill_from_root_file(bbfile, "tree", keys(bbfile["tree"]))              |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!

binWidth = 1
with(   
        title = "2D Vertex separation; "* L"d = \sqrt{(y_1 - y_2)^2 + (z_1 - z_2)^2}",
        legend = :best,
        thickness_scaling = 1.4,
        lw = 4,
        right_margin=12Plots.mm,
        bottom_margin=4Plots.mm,
        size = (1200, 800),
        nbins = 0:binWidth:100,
        normed = :true,

    ) do
    stephist(Bi214.d, c= 2,  label = "Bi214", lw = 3)
    stephist!(Tl208.d, c = 3,  label = "Tl208", lw = 3)
    stephist!(K40.d, c= 4, label = "K40", lw = 3)
    stephist!(Pa234m.d, c=5, label = "Pa234m", lw = 3)
    stephist!(
        bb.d, 
        c =1,
        label = "bb", 
        lw = 3, 
        fill = 0, 
        fa = 0.4, 
        fillstyle =:cross,
        ylims = (0, 0.04), 
        xlabel = "vertexSeparation/mm", 
        ylabel = "normalized rate/$binWidth mm"
    )
    safesave(plotsdir("VertexSeparation", "1Dhist_VertexSeparation_d.png"), current())
    current()
end

with(   
        title = "1D Vertex separation; "* L"d_y = |y_1 - y_2|",
        legend = :best,
        thickness_scaling = 1.4,
        lw = 4,
        right_margin=12Plots.mm,
        bottom_margin=4Plots.mm,
        size = (1200, 800),
        nbins = 0:binWidth:100,
        normed = :true,

    ) do
    stephist(Bi214.dy,  c= 2,label = "Bi214", lw = 3)
    stephist!(Tl208.dy,  c= 3, label = "Tl208", lw = 3)
    stephist!(K40.dy,  c=4, label = "K40", lw = 3)
    stephist!(Pa234m.dy, c=5 ,label = "Pa234m", lw = 3)
    stephist!(
        bb.dy, 
        c =1,
        label = "bb", 
        lw = 3, 
        fill = 0, 
        fa = 0.4, 
        ylims = (0, 0.08),
        fillstyle =:cross,
        xlabel = "vertexSeparation/mm", 
        ylabel = "normalized rate/$binWidth mm"
    )
    # safesave(plotsdir("VertexSeparation", "1Dhist_VertexSeparation_dy.png"), current())
    current()
end

with(   
        title = "1D Vertex separation; "* L"d_z = |z_1 - z_2|",
        legend = :best,
        thickness_scaling = 2.0,
        lw = 4,
        right_margin=12Plots.mm,
        bottom_margin=4Plots.mm,
        size = (1200, 800),
        nbins = 0:binWidth:100,
        normed = :true,
        frame = :true

    ) do
    stephist(Bi214.dz,  c= 2,label = "Bi214", lw = 3)
    stephist!(Tl208.dz,  c= 3, label = "Tl208", lw = 3)
    stephist!(K40.dz,  c=4, label = "K40", lw = 3)
    stephist!(Pa234m.dz, c=5,label = "Pa234m", lw = 3)
    stephist!(
        bb.dz, 
        c =1,
        label = "bb", 
        lw = 3, 
        fill = 0, 
        fa = 0.4, 
        ylims = (0, 0.05),
        fillstyle =:cross,
        xlabel = "vertexSeparation/mm", 
        ylabel = "normalized rate/$binWidth mm"
    )
    # safesave(plotsdir("VertexSeparation", "1Dhist_VertexSeparation_dz.png"), current())
    current()
end

quantile(bb.dy, 0.90)
quantile(bb.dz, 0.90)

with() do
    plot(
        Hist2D(
            (bb.dy, bb.dz),
            (0:1:50, 0:1:50),

        ),
        xlabel = "dy/mm",
        ylabel = "dz/mm",
        title = "2D histogram: bb",
        thickness_scaling = 1.9,
        right_margin = 10Plots.mm,
        size = (1000, 800)
    )
    # safesave(plotsdir("VertexSeparation", "2D_Histogram_bb.png"), current())
    current()
end

()
h(x, y) = exp(x^2 + y^2)

x = range(-3, 3, length=100)
y = range(-3, 3, length=100)
z = @. h(x', y)

tv = 0:8
tl = [L"10^{%$i}" for i in tv]
contourf(x, y, log10.(z), color=:turbo, levels=8,
    colorbar_ticks=(tv, tl), aspect_ratio=:equal,
    title=L"\exp(x^{2} + y^{2})", xlabel=L"x", ylabel=L"y")

pgfplotsx()