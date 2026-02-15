using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist
using UnROOT

println("loaded pkgs")

f = ROOTFile(datadir("data/final_phd/phase1_data.root"))
d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame

d.deltaCaloTime = abs.(d.caloTime1 .- d.caloTime2)

my_roi_1 = (
    sumE = (300, 2700),
    reconstructedEnergy1 = (350, 5000),
    reconstructedEnergy2 = (350, 5000),
    dy = (0, 100),
    dz = (0,100),
    deltaCaloTime = (0, 1.5)
)


function filter_data!(d, roi)
    d_new = deepcopy(d)
    for key in keys(roi)
        lo, hi = roi[key]
        filter!(row -> lo <= row[key] <= hi, d_new)
    end
    return d_new
end
d_roi_1 = filter_data!(d, my_roi_1)

h1 = Hist1D(d.sumE; binedges = (0:100:4000))
h2 = Hist1D(d_roi_1.sumE; binedges = (0:100:4000))

f = Figure(size = (1200, 800), fontsize = 30)
a1 = Axis(f[1,1], xlabel = "sumE")
a2 = Axis(f[2,1], xlabel = "sumE")

p1 = stephist!(a1, h1, label = "topology only")
p2 = stephist!(a2, h2, label = "sept cuts")

axislegend(a1)
axislegend(a2)
f

