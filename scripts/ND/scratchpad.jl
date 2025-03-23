using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, Combinatorics, LinearAlgebra
using Revise
using BlackBoxOptim

# Revise.track(SensitivityModule)

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(scriptsdir("Params.jl"))


vars = [
    # "phi", 
    "sumE", 
    # "maxE", 
    # "minE", 
    # "r", 
    "dy", 
    "dz"
    ]

bins = (
    # phi = (0,180),
    sumE = (0, 3500),
    # maxE = (0, 3000),
    # minE = (0, 3500),
    # r = (0, 100),
    dy = (-100, 100),
    dz = (-100, 100)
)


processes = load_ndim_processes("fal5_8perc_Boff_TIT", bins, vars)

signal = get_process("bb0nu_foil_bulk", processes)
# signal = get_process("bb0nuM1_foil_bulk", processes)
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", processes),
    get_process("Bi214_foil_bulk", processes),
    get_process("Bi214_wire_surface", processes),
    get_process("Tl208_foil_bulk", processes),
    # get_process("K40_foil_bulk", processes),
    # get_process("Pa234m_foil_bulk", processes),
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!( signal, 1e8 )
# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 0.99e8 )
set_nTotalSim!( background[2], 1e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 1e8 )
# set_nTotalSim!( background[5], 1e8 )
# set_nTotalSim!( background[6], 1e8 )

α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


prob(x) = - SensitivityModule.get_s_to_b(SNparams, α, vcat(signal, background), x;
    approximate="table")


function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


searchRange = make_stepRange(signal)
x0 = [
    # rand(range(bins.phi[1], bins.phi[2], 100)), rand(range(bins.phi[1], bins.phi[2], 100)), 
    rand(range(bins.sumE[1], bins.sumE[2], 100)), rand(range(bins.sumE[1], bins.sumE[2], 100)), 
    # rand(range(bins.maxE[1], bins.maxE[2], 100)), rand(range(bins.maxE[1], bins.maxE[2], 100)), 
    # rand(range(bins.minE[1], bins.minE[2], 100)), rand(range(bins.minE[1], bins.minE[2], 100)), 
    # rand(range(bins.r[1], bins.r[2], 100)), rand(range(bins.r[1], bins.r[2], 100)), 
    rand(range(bins.dy[1], bins.dy[2], 100)), rand(range(bins.dy[1], bins.dy[2], 100)), 
    rand(range(bins.dz[1], bins.dz[2], 100)), rand(range(bins.dz[1], bins.dz[2], 100)), 
    ]
res = bboptimize(
    prob,
    x0; 
    SearchRange = searchRange, 
    NumDimensions = length(searchRange),
    Method=:adaptive_de_rand_1_bin, 
    MaxTime = 60,#6*3600,
    TraceMode = :compact
)

# best_res = [8.17838, 180.0, 2702.3, 3104.02, 1352.34, 2974.99, 28.5889, 2218.03, 0.00107271, 82.4988, -81.4286, 97.2736, -96.9858, 68.2162]

# optim_res = Optim.optimize(prob, x0, NelderMead(), autodiff=:forward)
# optim_res.minimizer

function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    # best = res
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

best_roi =  get_best_ROI_ND(res, signal)

best_sens = get_sensitivityND(
    SNparams, 
    α, 
    vcat(signal, background), 
    best_roi
)

println(best_sens)

