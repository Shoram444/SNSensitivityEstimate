using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, LinearAlgebra
using BlackBoxOptim


# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(scriptsdir("Params.jl"))

# tt = collect(t)
# DataProcessND( tt,bins,vars,singleEParams[:bb_foil_bulk])

# mutable struct tester1
#     data::Vector{UnROOT.LazyEvent}
#     bins::NamedTuple
#     varNames::Vector{String}
#     singleEParams::Dict
# end

# tester1(tt,bins,vars,singleEParams[:bb_foil_bulk])

vars = [
    "phi", 
    "sumE", 
    "maxE", 
    # "minE", 
    "r", 
    # "dy", 
    # "dz"
    ]

bins = (
    phi = (0,180),
    sumE = (0, 3500),
    maxE = (0, 3500),
    # minE = (0, 3500),
    r = (0, 100),
    # dy = (-100, 100),
    # dz = (-100, 100)
)


processes = load_ndim_processes("fal5_TKrec", bins, vars)

# signal = get_process("bb0nu_foil_bulk", processes)
signal = get_process("bb0nuM1_foil_bulk", processes)
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", processes),
    get_process("Bi214_foil_bulk", processes),
    get_process("Bi214_wire_surface", processes),
    get_process("Tl208_foil_bulk", processes),
    get_process("K40_foil_bulk", processes),
    get_process("Pa234m_foil_bulk", processes),
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
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )

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
x0 = [
    rand(range(bins.phi[1], bins.phi[2], 100)), rand(range(bins.phi[1], bins.phi[2], 100)), 
    rand(range(bins.sumE[1], bins.sumE[2], 100)), rand(range(bins.sumE[1], bins.sumE[2], 100)), 
    rand(range(bins.maxE[1], bins.maxE[2], 100)), rand(range(bins.maxE[1], bins.maxE[2], 100)), 
    # rand(range(bins.minE[1], bins.minE[2], 100)), rand(range(bins.minE[1], bins.minE[2], 100)), 
    rand(range(bins.r[1], bins.r[2], 100)), rand(range(bins.r[1], bins.r[2], 100)), 
    # rand(range(bins.dy[1], 0.0, 100)), rand(range(0.0, bins.dy[2], 100)), 
    # rand(range(bins.dz[1], 0.0, 100)), rand(range(0.0, bins.dz[2], 100)), 
    ]


prob(x0)
prob(float.([0,10, 0,10, 0,10, 0,10]))
@time prob(float.([0,180, 0,3100, 0,2500, 0,100]))

res = bboptimize(
    prob,
    x0; 
    SearchRange = searchRange, 
    NumDimensions = length(searchRange),
    Method=:adaptive_de_rand_1_bin_radiuslimited, 
    MaxTime = 5*60,#24*3600,
    # TraceMode = :silent,
    MutationRate = 0.8, # Increase mutation to escape plateaus
    CrossoverRate = 0.9, # Allow more recombination of solutions
    PopulationSize = 200  # Increase population size for more diversity
)

function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
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

