using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, Combinatorics, LinearAlgebra
using Revise
using BlackBoxOptim

Revise.track(SensitivityModule)

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
include(scriptsdir("Params.jl"))


vars = ["phi", "reconstructedEnergy1"]

bins = (
    phi = (0,180),
    reconstructedEnergy1 = (0,1000),
)

function load_ndim_processes(dir::String, bins::NamedTuple, varNames::Vector{String})
    include(scriptsdir("Params.jl"))
    full_dir = datadir("sims", dir)
    processes = DataProcessND[]

    println("Loading files from: $full_dir ...")
    println("mode: NDim ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        if( split(file, ".")[end] != "root" )
            continue
        end
        
        f = ROOTFile(joinpath(full_dir, file)) 
        if(!haskey(f, "tree"))
            continue
        end

        data = LazyTree(f, "tree", varNames)

        fileName = split(file, ".")[1]  |> split |> first 
        push!(
                processes,
                DataProcessND(
                    data,
                    bins,
                    vars,
                    singleEParams[Symbol(fileName)]
                )
            )
        println("$fileName loaded")

    end
    println("Loaded $nFiles files.")

    return processes
end

processes = load_ndim_processes("fal5_8perc_Boff_TIT_twoDistinct_edep_bcu", bins, vars)

signal = get_process("bb0nu_foil_bulk", processes)
# signal = get_process("bb0nuM1_foil_bulk", data_processes)
# signal = get_process("bb0nuM2_foil_bulk", data_processes)

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
set_nTotalSim!( signal, 0.98e8 )
# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 0.99e8 )
set_nTotalSim!( background[2], 0.96e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 0.76e8 )
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )

@time get_roi_effciencyND(processes[1], bins).eff

@time get_roi_bkg_counts(
    processes, 
    bins
)


α= 1.64485362695147

@time get_s_to_b(
    SNparams, 
    α, 
    vcat(signal, background), 
    [0,180,0,1000];
    approximate="table"
)

prob(x) = - SensitivityModule.get_s_to_b(SNparams, α, vcat(signal, background), x)

function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end

searchRange = make_stepRange(signal)

res = bboptimize(
    prob; 
    SearchRange = searchRange, 
    NumDimensions = length(searchRange),
    Method=:adaptive_de_rand_1_bin, 
    MaxTime = 5*60,
    # InitialPopulation = [[0, 180, 400, 3000, 2700, 3200]]
)

function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

best_roi = get_best_ROI_ND(res, signal)

get_sensitivityND(
    SNparams, 
    α, 
    vcat(signal, background), 
    best_roi
)