using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
# using SNSensitivityEstimate, 
using CairoMakie, DataFramesMeta, CSV, Random, FHist
using UnROOT

println("loaded pkgs")

dur_1, dur_2, dur_3 = 25.97864583333333, 92.29909375, 52.5251875
tot_dur = dur_1 + dur_2 + dur_3

function read_data(file)
    f = ROOTFile(file)
    d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
    d.deltaCaloTime = abs.(d.caloTime1 .- d.caloTime2)
    return d
end

function get_n_exp(d; duration_seconds = tot_dur * 24*3600, nTotalSim = 1e8, activity = 1e-3)
    eff = nrow(d) / nTotalSim
    return eff * activity * duration_seconds
end


d1 = read_data(datadir("data/final_phd/phase1_data.root"))
d2 = read_data(datadir("data/final_phd/phase2_data.root"))
d3 = read_data(datadir("data/final_phd/phase3_data.root"))

d_sim_bb = read_data(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/bb_foil_bulk.root"))

free_roi = (
    sumE = (0, 5000),
    phi = (0,180)
)

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

function get_combined_data_hist(d1,d2,d3, roi, dur1, dur2, dur3)
    d1_roi = filter_data!(d1, roi)
    d2_roi = filter_data!(d2, roi)
    d3_roi = filter_data!(d3, roi)

    h1 = Hist1D(d1_roi.sumE; binedges = (0:100:4000)) * inv(dur1)
    h2 = Hist1D(d2_roi.sumE; binedges = (0:100:4000)) * inv(dur2)
    h3 = Hist1D(d3_roi.sumE; binedges = (0:100:4000)) * inv(dur3)
    return h1, h2, h3
end

h1, h2, h3 = get_combined_data_hist(d1, d2, d3, free_roi, dur_1, dur_2, dur_3)
h1_cut, h2_cut, h3_cut = get_combined_data_hist(d1, d2, d3, my_roi_1, dur_1, dur_2, dur_3)

d_sim_bb_roi = filter_data!(d_sim_bb, my_roi_1)
n_bb_exp = get_n_exp(d_sim_bb_roi, duration_seconds = tot_dur * 24*3600, nTotalSim = 1e8, activity = 1e-2) / tot_dur
h_sim_bb = normalize(Hist1D(d_sim_bb_roi.sumE; binedges = (0:100:4000)), width = false) * n_bb_exp 

f = Figure(size = (1200, 800), fontsize = 30)
a1 = Axis(f[1,1], xlabel = "sumE")
a2 = Axis(f[2,1], xlabel = "sumE")

p1 = stephist!(a1, h1+h2+h3, label = "topology only")
p2 = stephist!(a2, h1_cut+h2_cut+h3_cut, label = "sept cuts")
p3 = stephist!(a2, h_sim_bb, label = "bb foil bulk")

axislegend(a1)
axislegend(a2)
f
save("phase1_test_cuts.png", f)

