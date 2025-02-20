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


mutable struct NDimROIEfficiency
    binESingle
    binESum
    binPhi
    eff
end

import Base: *

*(x::Float64, e::NDimROIEfficiency) = x * e.eff

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
    range(0, 3500, 8),
    range(0, 3500, 8),
    1.0
)

f2 = ROOTFile("data/sims/fal5_8perc_Boff_TIT_twoDistinct_edepbcu/K40_foil_bulk.root")
d2 = LazyTree(f2, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])
phi2 = collect(d2.phi)
e2 = [maximum([ev.reconstructedEnergy1, ev.reconstructedEnergy2]) for ev in d2]
es2 = [sum([ev.reconstructedEnergy1, ev.reconstructedEnergy2]) for ev in d2]

p2 = NDimDataProcess(
    phi2,
    e2,
    es2,
    "K40",
    false,
    1.0,
    1.0,
    1e8,
    range(0, 180, 2),
    range(0, 3500, 8),
    range(0, 3500, 8),
    1.0
)

function get_bkg_counts(processes::Vector{NDimDataProcess})
    bkg = zeros(length(processes[1].eff))
    
    for p in processes
        p.signal && continue
        bkg += [p.amount * e.eff * p.activity * p.timeMeas for e in p.eff]
    end
    return bkg
end


function get_sensitivities(SNparams, α, processes::Vector{NDimDataProcess}; approximate="formula")
    signal_id = findall([p.signal for p in processes])
    length(signal_id) > 1 && @error "Only one signal process allowed! Provided $(length(signal_id))"
    
    signal_process = processes[first(signal_id)]
    ε = signal_process.eff
    b = get_bkg_counts(processes)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
    S_b = get_FC.(b, α; approximate=approximate)

    @show tHalf = @. constantTerm * ε / S_b
    length(tHalf)
    sensitivities = Vector{NDimSensitivityEstimate}(undef, length(tHalf))
    for i in 1:length(tHalf)
        sensitivities[i] = NDimSensitivityEstimate(
            signal_process.eff[i].binESingle,
            signal_process.eff[i].binESum,
            signal_process.eff[i].binPhi,
            tHalf[i],
            ε[i].eff,
            b[i]
        )
    end
    return sensitivities
end

α = 1.65
sensitivities = get_sensitivities(SNparams, α, [p1,p2]; approximate="formula")

sensitivities[1]

import Base: maximum
maximum(sensitivities::Vector{NDimSensitivityEstimate}) = sensitivities[argmax([s.tHalf for s in sensitivities])]

maximum(sensitivities)
