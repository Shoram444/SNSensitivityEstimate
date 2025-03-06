using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes, SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, BlackBoxOptim, Optim
using Revise
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

#mock signal, background data
analysisDict = Dict(
    :Bfield => "Boff", # magnetic field on/off
    :Eres => "8perc", # FWHM of the OMs (sorry for the naming...)
    :mode => "sumE", 
    :trackAlgo => "TIT",
    :signal => "bb0nu",
    :neutron_config => "full_shielding"
)

files_directory = "fal5_$(analysisDict[:Eres])_$(analysisDict[:Bfield])_$(analysisDict[:trackAlgo])_twoDistinct_edep_bcu"

binsAngle = 0:5:180
binsESingle = 0:100:3500
binsESum = 0:100:3500
# Load all the processes in the directory. Function `load_processes` takes two arguments:
#  1. dir::String -> the name of the directory where the root files are stored
#  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
data_processes = load_ndim_processes(
    files_directory, binsAngle, binsESingle, binsESum
)



signal = get_process("bb0nuM2_foil_bulk", data_processes)

background = [
    get_process("bb_foil_bulk", data_processes),
    get_process("Bi214_foil_bulk", data_processes),
    get_process("Bi214_wire_surface", data_processes),
    get_process("Tl208_foil_bulk", data_processes),
    get_process("K40_foil_bulk", data_processes),
    get_process("Pa234m_foil_bulk", data_processes),
    # get_process("neutron_external", hist_processes, "current shielding")
]

# set_nTotalSim!(signal, 0.98e8)
set_nTotalSim!(signal, 1e8)

set_nTotalSim!(background[1], 0.99e8)
set_nTotalSim!(background[2], 0.96e8)
set_nTotalSim!(background[3], 1e8)
set_nTotalSim!(background[4], 0.76e8)
set_nTotalSim!(background[5], 1e8)
set_nTotalSim!(background[6], 1e8)

set_signal!(background[1], false)

α = 1.64485362695147

passes_roi(signal.data, signal.bins, signal.bins)

f3d(x) = -get_s_to_b(SNparams, α, vcat(signal, background), x)

function make_stepRange(process)
    a = (process.binsAngle[1], process.binsAngle[end])
    es = (process.binsESingle[1], process.binsESingle[end])
    eS = (process.binsESum[1], process.binsESum[end])
    return [a, a, es, es, eS, eS]
end
searchRange = make_stepRange(signal)

res = bboptimize(
    f3d; 
    SearchRange = searchRange, 
    NumDimensions = 6,
    Method=:adaptive_de_rand_1_bin, 
    MaxTime = 1*60,
    InitialPopulation = [[0, 180, 400, 3000, 2700, 3200]]
)

function get_best_ROI3D(res)
    best = best_candidate(res)
    best_roi = Dict(
        :angle => round.(best[1:2], digits =-1),
        :esingle => round.(best[3:4], digits =-2),
        :esum => round.(best[5:6], digits =-2)
    )
    return best_roi
end

best_roi = get_best_ROI3D(res)
@profview best_sens= get_sensitivity3D(
    SNparams, 
    α, 
    vcat(signal, background), 
    best_roi;
    approximate="table"
)


let
    f = Figure()
    axAngle = Axis(f[1,1], xlabel = L"angle ($^{\circ}$)")
    p1 = stephist!(axAngle, ps.dataAngle, bins=180, label="0nu", normalization = :pdf)
    stephist!(axAngle, 180 .- ps.dataAngle, bins=180, label="neutron", normalization = :pdf)
    stephist!(axAngle, pb.dataAngle, bins=180, label="bkg", normalization = :pdf)
    axislegend(axAngle)
    # axEsingle = Axis(f[1,2])
    # p2 = stephist!(axEsingle, ps.dataESingle, bins=ps.binsESingle, label="signal", normalization = :pdf)
    # stephist!(axEsingle, pb.dataESingle, bins=ps.binsESingle, label="signal", normalization = :pdf)
    f
end



let 
    f = Figure()
    ax = Axis(f[1,1], xlabel = "b", ylabel = "S(b)")
    p1 = lines!(ax, 0:1:350, map(x -> 1.8 * sqrt(x), 0:1:350), label = "1.8 * √b")
    p1 = lines!(ax, 0:1:350, map(x -> α * sqrt(x), 0:1:350), label = "α * √b")
    p1 = lines!(ax, 0:1:350, map(x -> get_FC(x, α; approximate="formula"), 0:1:350), label = "FC approx. formula")
    p1 = lines!(ax, 0:1:350, map(x -> get_FC(x, α; approximate="table"), 0:1:350), label = "FC lookup table")
    axislegend(ax, position = :rb)
    save("FC_comparison.png", f)
    f
end


