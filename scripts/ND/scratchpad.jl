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

mutable struct SensitivityEstimateND
    roi::NamedTuple
    tHalf
    signalEff
    bkgCounts
end

vars = ["phi", "reconstructedEnergy1"]

bins = (
    phi = (0,180),
    reconstructedEnergy1 = (0,1000),
)


mutable struct ROIEfficiencyND
    roi::NamedTuple
    eff
end


import Base: *, maximum
*(x::Float64, e::ROIEfficiencyND) = x * e.eff

mutable struct DataProcessND <: AbstractProcess
    data::LazyTree
    isotopeName::String
    signal::Bool
    activity::Real
    timeMeas::Real
    nTotalSim::Real
    bins::NamedTuple
    amount::Real
end

using Polyester

function get_roi_effciencyND(
    data::LazyTree, 
    roi::NamedTuple, 
    nTotalSim::Real
)
    k_ = keys(roi)  # Extract keys once to avoid repeated lookup

    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@thread for i in eachindex(data)
        if passes_roi( data[i], roi, k_)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end

    return ROIEfficiencyND(
        roi,  # Assuming bins should be roi (you can replace this if incorrect)
        count[] / nTotalSim  # Extract final atomic count
    )
end

f = ROOTFile(datadir("sims", "fal5_8perc_Boff_TIT_twoDistinct_edep_bcu", "Bi214_foil_bulk.root"))
lt1 = LazyTree(f, "tree", vars)

f2 = ROOTFile(datadir("sims", "fal5_8perc_Boff_TIT_twoDistinct_edep_bcu", "bb_foil_bulk.root"))
lt2 = LazyTree(f2, "tree", vars)

@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::NamedTuple,
    k_
)
    @inbounds for k in k_
        val = data[k]  # Avoid repeated dictionary lookups
        range = roi[k]  # Avoid repeated access to NamedTuple
        if !(range[1] ≤ val < range[2])
            return false
        end
    end
    return true
end

function get_roi_effciencyND(
    process::DataProcessND,
    roi::NamedTuple
)
    return get_roi_effciencyND(
        process.data,
        roi,
        process.nTotalSim
    )
end

test_bins = (phi=(0,180), reconstructedEnergy1=(0,1000))
@time passes_roi(lt1[1], test_bins, keys(test_bins))
@time get_roi_effciencyND(lt1, test_bins, 1e8)

function DataProcessND(
    data::LazyTree,
    bins::NamedTuple,
    processDict::Dict
)
    @unpack isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount = processDict
    println("creating process: $isotopeName")

    return DataProcessND(
        data,
        isotopeName,
        signal,
        activity,
        timeMeas,
        nTotalSim,
        bins,
        amount
    )
end

ps = DataProcessND(
    lt2,
    bins,
    Dict(
        :isotopeName => "test_process",
        :signal => true,
        :activity => 1.0,
        :timeMeas => 1.0,
        :nTotalSim => 1.0,
        :bins => bins,
        :amount => 1.0
    )
)

pb = DataProcessND(
    lt1,
    bins,
    Dict(
        :isotopeName => "bkg",
        :signal => false,
        :activity => 1.0,
        :timeMeas => 1.0,
        :nTotalSim => 1.0,
        :bins => bins,
        :amount => 1.0
    )
)

function get_roi_bkg_counts(
    processes::Vector{DataProcessND}, 
    roi::NamedTuple
)
    bkg = 0.0

    for p in processes
        p.signal && continue
        eff = get_roi_effciencyND(p, roi).eff
        bkg += p.amount * eff * p.activity * p.timeMeas
    end
    return bkg
end

@time get_roi_effciencyND(ps, test_bins).eff

@time get_roi_bkg_counts(
    [ps, pb], 
    test_bins
)


function get_s_to_b1(SNparams, 
    α, 
    processes::Vector{DataProcessND}, 
    ROIs::Vector{<:Real};
    approximate="table"
)

    bin_keys = keys(processes[1].bins)

    # Construct ROI 
    roi = NamedTuple(k => (round(ROIs[i]), round(ROIs[i+1])) 
                     for (i, k) in zip(1:2:length(bin_keys)*2-1, bin_keys))

    signal_id = findfirst(p -> p.signal, processes)
    if signal_id === nothing
        @error "No signal process found!"
        return 0.0
    end
                    
    signal_process = processes[signal_id]
    ε = get_roi_effciencyND(signal_process.data, roi, signal_process.nTotalSim).eff
    ε == 0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return ε / S_b
end


function get_sensitivityND(
    SNparams, 
    α, 
    processes::Vector{DataProcessND}, 
    roi::NamedTuple;
    approximate="table"
)

    #check that rois are within the range of the data
    for k in keys(roi)
        if (roi[k][1] < processes[1].bins[k][1] || roi[k][2] > processes[1].bins[k][2])
            @error "ROI out of range for $k"
        end
    end

    signal_id = findall([p.signal for p in processes])
    length(signal_id) > 1 && @error "Only one signal process allowed! Provided $(length(signal_id))"

    signal_process = processes[first(signal_id)]
    ε = get_roi_effciencyND(signal_process.data, roi, signal_process.nTotalSim).eff
    ε == 0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
    S_b = get_FC(b, α; approximate=approximate)

    tHalf = constantTerm * ε / S_b

    # return tHalf
    return SensitivityEstimateND(
            roi,
            tHalf,
            ε,
            b
        )

end

α= 1.64485362695147
get_sensitivityND(
    SNparams, 
    α, 
    [ps, pb], 
    (phi=(0,180), reconstructedEnergy1=(0,1000));
    approximate="table"
)

prob(x) = -get_s_to_b1(SNparams, α, vcat(ps, pb), x)

function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end

searchRange = make_stepRange(ps)

res = bboptimize(
    prob; 
    SearchRange = searchRange, 
    NumDimensions = length(searchRange),
    Method=:adaptive_de_rand_1_bin, 
    MaxTime = 1*60,
    # InitialPopulation = [[0, 180, 400, 3000, 2700, 3200]]
)

function get_best_ROI3D(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

best_roi = get_best_ROI3D(res, ps)
