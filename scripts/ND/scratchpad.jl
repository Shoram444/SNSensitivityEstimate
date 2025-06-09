using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, LinearAlgebra
using BlackBoxOptim
Revise.track(SensitivityModule)

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(srcdir("params/Params.jl"))

vars = [
    "phi", 
    "sumE", 
    # "maxE", 
    # "minE", 
    "r", 
    # "dy", 
    # "dz",
    # "sameSide"
    "lPint",
    "lPext"
    ]

bins = (
    phi = (0,180),
    sumE = (300, 3500),
    # maxE = (0, 3500),
    # minE = (0, 3500),
    r = (0, 200),
    # dy = (0, 300),
    # dz = (0, 300),
    # sameSide = (0, 1)
    lPint = (0, 100),
    lPext = (0, 100)
)

processes = load_ndim_processes("fal5_TKrec_J40", bins, vars)

signal = get_process("bb0nu_foil_bulk", processes) |> first
# signal = get_process("RH037_foil_bulk", processes) |> first
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", processes) |> first,
    get_process("Bi214_foil_bulk", processes) |> first,
    get_process("Bi214_wire_surface", processes) |> first,
    get_process("Tl208_foil_bulk", processes) |> first,
    get_process("K40_foil_bulk", processes) |> first,
    get_process("Pa234m_foil_bulk", processes) |> first,
    get_process("gamma_experimental_surface", processes) |> first,
    # get_process("Bi214_PMT_glass_bulk", processes) |> first,
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!( signal, 0.1e8 )
# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 0.1e8 )
set_nTotalSim!( background[2], 1e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 1e8 )
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )
set_nTotalSim!( background[7], 5e8 )

α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


prob(x) = - SensitivityModule.get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="formula")


function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


searchRange = make_stepRange(signal)
# x0 = [
#     # rand(range(bins.phi[1], bins.phi[2], 100)), rand(range(bins.phi[1], bins.phi[2], 100)), 
#     # rand(range(bins.sumE[1], bins.sumE[2], 100)), rand(range(bins.sumE[1], bins.sumE[2], 100)), 
#     # rand(range(bins.maxE[1], bins.maxE[2], 100)), rand(range(bins.maxE[1], bins.maxE[2], 100)), 
#     # rand(range(bins.minE[1], bins.minE[2], 100)), rand(range(bins.minE[1], bins.minE[2], 100)), 
#     rand(range(bins.r[1], bins.r[2], 100)), rand(range(bins.r[1], bins.r[2], 100)), 
#     # rand(range(bins.dy[1], 0.0, 100)), rand(range(0.0, bins.dy[2], 100)), 
#     # rand(range(bins.dz[1], 0.0, 100)), rand(range(0.0, bins.dz[2], 100)), 
#     ]

x0 = float.([
        0, 
        180, 
        2700, 
        3100, 
        0, 
        50, 
        0, # pint max 
        1.4, # pint min
        3,  # pext min
        100 # pext max
    ])

using Optim

lower = [0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0]
upper = [180.0, 180.0, 3500.0, 3500.0, 50.0, 50.0, 100.0, 100.0, 100.0, 100.0]

# result = optimize(prob, lower, upper, x0, Fminbox(BFGS()), Optim.Options(show_trace = true, iterations = 5))
# best2 = result.minimizer
# best2[2] = 180.0
# prob(x0)

# prob(float.([0,10, 0,10, 0,10, 0,10]))
# @time prob(float.([0,180, 0,3100, 0,2500, 0,100]))


res = bboptimize(
    prob,
    x0; 
    SearchRange = searchRange, 
    NumDimensions = length(searchRange),
    Method=:adaptive_de_rand_1_bin_radiuslimited, 
    MaxTime = 0.5*60,#24*3600,
    TraceMode = :compact,
    PopulationSize = 800  # Increase population size for more diversity
)

function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => ((best[i]), (best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
    best_roi = NamedTuple(
        k => ((best[i]), (best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

@show best_roi =  get_best_ROI_ND(res, signal)

# best_sens = get_sensitivityND(
#     SNparams, 
#     α, 
#     vcat(signal, background), 
#     best_roi
# ) |> println
best_sens = get_sensitivityND(
    SNparams, 
    α, 
    vcat(signal, background), 
    best_roi
) |> println


let
    res2 = float.([
        0, 
        180, 
        2700, 
        3100, 
        0, 
        50, 
        0, # pint max 
        1.4, # pint min
        3,  # pext min
        100 # pext max
        ])
    best2 = get_best_ROI_ND(res2, signal)
    println(get_sensitivityND(SNparams, α, vcat(signal, background), best2; approximate="table", add_mock_bkg=0.0).tHalf)

end



# 3.46e24 with NEMO3 tof for 0nu