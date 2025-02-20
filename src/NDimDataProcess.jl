import Base: *, maximum

"""
    aa
"""
mutable struct NDimSensitivityEstimate
    binESingle
    binESum
    binPhi
    tHalf
    signalEff
    bkgCounts
end


"""
    Holds information about ROI efficiency
"""
mutable struct NDimROIEfficiency
    binESingle
    binESum
    binPhi
    eff
end


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

function NDimDataProcess(
        dataAngle::Vector{<:Real},
        dataESingle::Vector{<:Real},
        dataESum::Vector{<:Real},
        binsAngle,
        binsESingle,
        binsESum,
        processDict::Dict
    )
    @unpack isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount = processDict

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

    tHalf = @. constantTerm * ε / S_b
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


maximum(sensitivities::Vector{NDimSensitivityEstimate}) = sensitivities[argmax([s.tHalf for s in sensitivities])]
