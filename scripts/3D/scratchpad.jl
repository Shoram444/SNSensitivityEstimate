using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, Combinatorics, LinearAlgebra
using Revise
Revise.track(SensitivityModule)

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
include(scriptsdir("Params.jl"))


mutable struct NDimSensitivityEstimate
    binESingle
    binESum
    binPhi
    tHalf
    signalEff
    bkgCounts
end

function NDimSensitivityEstimate(binESingle, binESum, binPhi, tHalf, signalEff, bkgCounts)
    NDimSensitivityEstimate(binESingle, binESum, binPhi, tHalf, signalEff, bkgCounts)
end

mutable struct NDimROIEfficiency
    binESingle
    binESum
    binPhi
    eff
end

mutable struct NDimDataProcess <: AbstractProcess
    dataAngle::Vector{<:Real}
    dataESingle::Vector{<:Real}
    dataESum::Vector{<:Real}
    isotopeName::String
    signal::Bool
    activity::Real
    timeMeas::Real
    nTotalSim::Real
    binsAngle::AbstractRange
    binsESingle::AbstractRange
    binsESum::AbstractRange
    amount::Real
    eff::Vector{NDimROIEfficiency}
end

function make_ROI_combinations(bins::AbstractRange)
    return collect(combinations(bins, 2))
end


function get_nDim_effciency(
        dataAngle::Vector{<:Real}, 
        dataESingle::Vector{<:Real}, 
        dataESum::Vector{<:Real}, 
        binsAngle::Vector{<:Real}, 
        binsESingle::Vector{<:Real}, 
        binsESum::Vector{<:Real}, 
        nTotalSim::Real
    )
    nPassed = 0
    for i in 1:length(dataAngle)
        if (
            binsAngle[1] <= dataAngle[i] < binsAngle[2] && 
            binsESingle[1] <= dataESingle[i] < binsESingle[2] && 
            binsESum[1] <= dataESum[i] < binsESum[2]
        )
            nPassed += 1
        end
    end
    return nPassed/nTotalSim
end


function NDimDataProcess(
        dataAngle::Vector{<:Real},
        dataESingle::Vector{<:Real},
        dataESum::Vector{<:Real},
        isotopeName::String,
        signal::Bool,
        activity::Real,
        timeMeas::Real,
        nTotalSim::Real,
        binsAngle::AbstractRange,
        binsESingle::AbstractRange,
        binsESum::AbstractRange,
        amount::Real,
    )
    comb_angle = make_ROI_combinations(binsAngle)
    comb_esingle = make_ROI_combinations(binsESingle)
    comb_esum = make_ROI_combinations(binsESum)

    all_combinations = collect(Iterators.product(comb_angle, comb_esingle, comb_esum))

    eff = Vector{NDimROIEfficiency}(undef, length(all_combinations))
    for (i, comb) in enumerate(all_combinations)
        ROIAngle = comb[1]
        ROIESingle = comb[2]
        ROIESum = comb[3]
         
        _eff = NDimROIEfficiency(
            ROIAngle,
            ROIESingle,
            ROIESum,
            get_nDim_effciency(dataAngle, dataESingle, dataESum, ROIAngle, ROIESingle, ROIESum, nTotalSim)
        )
        eff[i]  = _eff
    end

    return NDimDataProcess(dataAngle, dataESingle, dataESum, isotopeName, signal, activity, timeMeas, nTotalSim, binsAngle, binsESingle, binsESum, amount, eff)
end


f1 = ROOTFile("data/sims/fal5_8perc_Boff_TIT_twoDistinct_edepbcu/bb_foil_bulk.root")
d1 = LazyTree(f1, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])
phi1 = collect(d1.phi)
e1 = [maximum([ev.reconstructedEnergy1, ev.reconstructedEnergy2]) for ev in d1]
es1 = [sum([ev.reconstructedEnergy1, ev.reconstructedEnergy2]) for ev in d1]

p1 = NDimDataProcess(
    phi1,
    e1,
    es1,
    "bb",
    true,
    1.0,
    1.0,
    1e8,
    range(0, 180, 2),
    range(500, 3500, 3),
    range(500, 3500, 3),
    1.0
)

f2 = ROOTFile("data/sims/fal5_8perc_Boff_TIT_twoDistinct_edepbcu/K40_foil_bulk.root")
d2 = LazyTree(f2, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])
phi2 = collect(d2.phi)
e2 = [maximum([ev.reconstructedEnergy1, ev.reconstructedEnergy2]) for ev in d2]
es2 = [sum([ev.reconstructedEnergy1, ev.reconstructedEnergy2]) for ev in d2]

p1 = NDimDataProcess(
    phi2,
    e2,
    es2,
    "K40",
    true,
    1.0,
    1.0,
    1e8,
    range(0, 180, 2),
    range(500, 3500, 3),
    range(500, 3500, 3),
    1.0
)



