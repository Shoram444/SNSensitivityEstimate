using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

Bi214File = ROOTFile("data/sims/Bi214_EnePhiDistPos_J23.root")
Bi214SurfFile = ROOTFile("data/sims/Bi214_wire_surface_EnePhiDistPos_J23.root")
Bi214BulkFile = ROOTFile("data/sims/Bi214_wire_bulk_EnePhiDistPos_J23.root")
Tl208File = ROOTFile("data/sims/Tl208_EnePhiDistPos_J23.root")
K40File = ROOTFile("data/sims/K40_EnePhiDistPos_J23.root")
Pa234mFile = ROOTFile("data/sims/Pa234m_EnePhiDistPos_J23.root")
bbFile = ROOTFile("data/sims/2nubb_EnePhiDistPos_J23.root")   # Standard 2nubb


Bi214 = fill_from_root_file(Bi214File, "tree", keys(Bi214File["tree"]))     |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
Bi214Surf = fill_from_root_file(Bi214SurfFile, "tree", keys(Bi214File["tree"]))     |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
Bi214Bulk = fill_from_root_file(Bi214BulkFile, "tree", keys(Bi214File["tree"]))     |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
Tl208 = fill_from_root_file(Tl208File, "tree", keys(Tl208File["tree"]))     |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
K40 = fill_from_root_file(K40File, "tree", keys(K40File["tree"]))           |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
Pa234m = fill_from_root_file(Pa234mFile, "tree", keys(Pa234mFile["tree"]))  |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!
bb = fill_from_root_file(bbFile, "tree", keys(bbFile["tree"]))              |> add_vertex_2D_separation_column! |> add_vertex_dy_separation_column! |> add_vertex_dz_separation_column!

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
    stephist(Bi214.d, c= 2,  label = "Bi214", lw = 4)
    # stephist!(Bi214Surf.d, c= 6,  label = "Bi214Surf", lw = 4)
    # stephist!(Bi214Bulk.d, c= 8,  label = "Bi214Bulk", lw = 4)
    stephist!(Tl208.d, c = 3,  label = "Tl208", lw = 4)
    stephist!(K40.d, c= 4, label = "K40", lw = 4)
    stephist!(Pa234m.d, c=5, label = "Pa234m", lw = 4)
    stephist!(
        bb.d, 
        c =1,
        label = "bb", 
        lw = 4, 
        fill = 0, 
        fa = 0.4, 
        fillstyle =:cross,
        ylims = (0, 0.04), 
        xlabel = "vertexSeparation/mm", 
        ylabel = "normalized rate/$binWidth mm",
        widen = :false,
        framestyle=:box,
        minorgrid = :true,
        thickness_scaling = 1.5
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
    # stephist!(Bi214Surf.dy,  c= 6,label = "Bi214Surf", lw = 3)
    # stephist!(Bi214Bulk.dy,  c= 8,label = "Bi214Bulk", lw = 3)
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
    safesave(plotsdir("VertexSeparation", "1Dhist_VertexSeparation_dy.png"), current())
    current()
end

with(   
        title = "1D Vertex separation; "* L"d_z = |z_1 - z_2|",
        legend = :best,
        thickness_scaling = 1.4,
        lw = 4,
        right_margin=12Plots.mm,
        bottom_margin=4Plots.mm,
        size = (1200, 800),
        nbins = 0:binWidth:100,
        normed = :true,
    ) do
    stephist(Bi214.dz,  c= 2,label = "Bi214", lw = 3)
    # stephist!(Bi214Surf.dz,  c= 6,label = "Bi214Surf", lw = 3)
    # stephist!(Bi214Bulk.dz,  c= 8,label = "Bi214Bulk", lw = 3)
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
        fillstyle = :cross,
        ylims = (0, 0.05),
        xlabel = "vertexSeparation/mm", 
        ylabel = "normalized rate/$binWidth mm"
    )
    safesave(plotsdir("VertexSeparation", "1Dhist_VertexSeparation_dz.png"), current())
    current()
end

quantile(Bi214Surf.dy, 0.90)
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

let 
    vertexDistances = range(0, 300, 100)
    sigFrac = Float64[]
    bkgFoilFrac = Float64[]
    bkgWireFrac = Float64[]

    bkgFoilDistances = vcat( 
        Bi214.d,
        Tl208.d,
        K40.d,
        Pa234m.d,
    )

    bkgWireDistances = vcat( 
        Bi214Surf.d,
        Bi214Bulk.d,
    )

    for d in vertexDistances
        sigPass = filter(x -> x .< d, bb.d)
        push!( sigFrac, length(sigPass)/length(bb.d) )

        bkgFoilPass = filter(x -> x .< d, bkgFoilDistances)
        push!( bkgFoilFrac, length(bkgFoilPass)/length(bkgFoilDistances) )

        bkgWirePass = filter(x -> x .< d, bkgWireDistances)
        push!( bkgWireFrac, length(bkgWirePass)/length(bkgWireDistances) )
    end


    plot( 
        vertexDistances, 
        sigFrac .* 100, 
        ylims = (0,100), 
        label = "signal",
        lw = 3, 
        widen = :false,
        framestyle=:box
    )
    plot!( 
        vertexDistances, 
        bkgFoilFrac .*100, 
        label = "background sources from foil", 
        lw = 3, 

    )

    plot!( 
        vertexDistances, 
        bkgWireFrac .*100, 
        label = "background sources from wires", 
        xlabel = "vertex distance cut-off [mm]",
        ylabel = "fraction passed [%]" ,
        title = "fraction of events that pass 2D Vertex cut",
        right_margin = 8Plots.mm,
        legend= :bottomright,
        lw = 3, 

    )
    plot!(
        vertexDistances,
        sigFrac ./ ( bkgFoilFrac .+ bkgWireFrac ) .* 100,
        label = "signal/background",
        lw = 3
    )

    safesave( plotsdir("VertexSeparation", "Fraction_pass_2DVertex_cut.png"), current())
    current()
end


let 
    vertexDistances = range(0, 300, 100)
    sigFrac = Float64[]
    bkgFoilFrac = Float64[]
    bkgWireFrac = Float64[]

    bkgFoilDistances = vcat( 
        Bi214.dy,
        Tl208.dy,
        K40.dy,
        Pa234m.dy,
    )

    bkgWireDistances = vcat( 
        Bi214Surf.dy,
        Bi214Bulk.dy,
    )

    for d in vertexDistances
        sigPass = filter(x -> x .< d, bb.dy)
        push!( sigFrac, length(sigPass)/length(bb.dy) )

        bkgFoilPass = filter(x -> x .< d, bkgFoilDistances)
        push!( bkgFoilFrac, length(bkgFoilPass)/length(bkgFoilDistances) )

        bkgWirePass = filter(x -> x .< d, bkgWireDistances)
        push!( bkgWireFrac, length(bkgWirePass)/length(bkgWireDistances) )
    end


    plot( 
        vertexDistances, 
        sigFrac .* 100, 
        ylims = (0,100), 
        label = "signal",
        framestyle = :box,
        lw = 3
    )
    plot!( 
        vertexDistances, 
        bkgFoilFrac .*100, 
        label = "background sources from foil", 
        lw = 3
    )

    plot!( 
        vertexDistances, 
        bkgWireFrac .*100, 
        label = "background sources from wires", 
        xlabel = "vertex distance cut-off [mm]",
        ylabel = "fraction passed [%]" ,
        title = "fraction of events that pass dy Vertex cut",
        right_margin = 8Plots.mm,
        legend= :bottomright,
        lw = 3

    )
    plot!(
        vertexDistances,
        sigFrac ./ ( bkgFoilFrac .+ bkgWireFrac ) .* 100,
        label = "signal/background",
        lw = 3
    )
    safesave( plotsdir("VertexSeparation", "Fraction_pass_dy_vertex_cut.png"), current())
    current()
end

let 
    vertexDistances = range(0, 300, 100)
    sigFrac = Float64[]
    bkgFoilFrac = Float64[]
    bkgWireFrac = Float64[]

    bkgFoilDistances = vcat( 
        Bi214.dz,
        Tl208.dz,
        K40.dz,
        Pa234m.dz,
    )

    bkgWireDistances = vcat( 
        Bi214Surf.dz,
        Bi214Bulk.dz,
    )

    for d in vertexDistances
        sigPass = filter(x -> x .< d, bb.dz)
        push!( sigFrac, length(sigPass)/length(bb.dz) )

        bkgFoilPass = filter(x -> x .< d, bkgFoilDistances)
        push!( bkgFoilFrac, length(bkgFoilPass)/length(bkgFoilDistances) )

        bkgWirePass = filter(x -> x .< d, bkgWireDistances)
        push!( bkgWireFrac, length(bkgWirePass)/length(bkgWireDistances) )
    end


    plot( 
        vertexDistances, 
        sigFrac .* 100, 
        ylims = (0,100), 
        label = "signal",
        lw = 3,
        framestyle =:box,
        minorgrid = :true,
        widen = :false
    )
    plot!( 
        vertexDistances, 
        bkgFoilFrac .*100, 
        label = "background sources from foil", 
        lw = 3
    )

    plot!( 
        vertexDistances, 
        bkgWireFrac .*100, 
        label = "background sources from wires", 
        xlabel = "vertex distance cut-off [mm]",
        ylabel = "fraction passed [%]" ,
        title = "fraction of events that pass dz Vertex cut",
        right_margin = 8Plots.mm,
        legend= :bottomright,
        lw = 3
    )
    plot!(
        vertexDistances,
        sigFrac ./ ( bkgFoilFrac .+ bkgWireFrac ) .* 100,
        label = "signal/background",
        lw = 3
    )
    safesave( plotsdir("VertexSeparation", "Fraction_pass_dz_vertex_cut.png"), current())
    current()
end


a = 2.0  # Semi-major axis
b = 1.0  # Semi-minor axis
θ = 0  # Angle of rotation (in radians)

# Function to parametrically represent the ellipse
function ellipse(t)
    x = a * cos(t)
    y = b * sin(t)
    return (x * cos(θ) - y * sin(θ), x * sin(θ) + y * cos(θ))
end

# Generate points on the ellipse
t = range(0, stop=2π, length=100)
points = [ellipse(angle) for angle in t]

# Plot the ellipse
plot(points, aspect_ratio=1, label="Ellipse")
title!("Ellipse in Julia")
xlabel!("X-axis")
ylabel!("Y-axis")
