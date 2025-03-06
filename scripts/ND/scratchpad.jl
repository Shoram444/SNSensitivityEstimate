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

bins = (
    angle = (0,180),
    esingle = (0,1000),
    esum = (0,1000)
)

se = SensitivityEstimateND(bins, 100, 1, 1) 

mutable struct ROIEfficiencyND
    roi::NamedTuple
    eff
end

ROIEfficiencyND(bins, 0.5)

import Base: *, maximum
*(x::Float64, e::ROIEfficiencyND) = x * e.eff

mutable struct DataProcessND <: AbstractProcess
    data::NamedTuple
    isotopeName::String
    signal::Bool
    activity::Real
    timeMeas::Real
    nTotalSim::Real
    bins::NamedTuple
    amount::Real
end

function get_roi_effciencyND(
    data::NamedTuple, 
    roi::NamedTuple, 
    nTotalSim::Real
)
    nthreads = Threads.nthreads()  # Get number of available threads
    local_counts = fill(0, nthreads)  # Create thread-local counters

    Threads.@threads for tid in 1:nthreads
        local_count = 0  # Local variable for counting
        @inbounds for i in tid:nthreads:length(data[1])  # Chunking work by thread ID
            if (
                passes_roi(
                    NamedTuple(k => data[k][i] for k in keys(data)),
                    roi
                )
            )
                local_count += 1  # Update local count
            end
        end
        local_counts[tid] = local_count  # Store count from this thread
    end

    return ROIEfficiencyND(
        bins,
        sum(local_counts) / nTotalSim
    )
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

@inline function passes_roi(
    data::NamedTuple,
    roi::NamedTuple, 
)
    for k in keys(roi)
        if !(roi[k][1] ≤ data[k] < roi[k][2])
            return false
        end
    end
    return true
end

passes_roi((angle=90, esingle=500, esum=500), (angle=(0,180), esingle=(0,1000), esum=(0,1000)))

function DataProcessND(
    data::NamedTuple,
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
    (angle=rand(0:180, 10), esingle=rand(0:1000, 10), esum=rand(0:1000, 10)),
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
    (angle=rand(0:180, 10), esingle=rand(0:1000, 10), esum=rand(0:1000, 10)),
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

get_roi_bkg_counts(
    [ps, pb], 
    (angle=(0,180), esingle=(0,1000), esum=(0,200))
)

function get_s_to_b(SNparams, 
    α, 
    processes::Vector{DataProcessND}, 
    ROIs::Vector{<:Real};
    approximate="table"
)

    # construct ROI from the ROIs vector (which was provided by gradient descent) and round to the 10s
    roi = NamedTuple(
            k => (round(ROIs[i]), round(ROIs[i+1])) 
            for (i,k) in zip(1:2:length(processes[1].bins)*2-1, keys(processes[1].bins))
        )

    signal_id = findall([p.signal for p in processes])
    length(signal_id) > 1 && @error "Only one signal process allowed! Provided $(length(signal_id))"
    
    signal_process = processes[first(signal_id)]
    ε = get_roi_effciencyND(signal_process.data, roi, signal_process.nTotalSim).eff
    ε == 0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return ε / S_b
end

get_s_to_b(SNparams, 
    1.64485362695147,
    [ps, pb], 
    [0, 180, 0, 1000, 0, 1000];
    approximate="table"
)

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
    (angle=(0,180), esingle=(0,1000), esum=(0,1000));
    approximate="table"
)

f(x) = -get_s_to_b(SNparams, α, vcat(ps, pb), x)

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
    f; 
    SearchRange = searchRange, 
    NumDimensions = length(searchRange),
    Method=:adaptive_de_rand_1_bin, 
    MaxTime = 3*60,
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