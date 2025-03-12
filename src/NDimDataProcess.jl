import Base: *, maximum

"""
    Holds information about sensitivity estimate
    Fields:
    - `roi`: NamedTuple of the ROI
    - `tHalf`: sesnitivity estimate
    - `signalEff`: Signal efficiency
    - `bkgCounts`: Background counts
"""
mutable struct SensitivityEstimateND
    roi::NamedTuple
    tHalf
    signalEff
    bkgCounts
end


"""
    Holds information about ROI efficiency
"""
mutable struct ROIEfficiencyND
    roi::NamedTuple
    eff
end


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
    varNames::Tuple
end



function get_roi_effciencyND(
    data::LazyTree, 
    roi::NamedTuple, 
    nTotalSim::Real,
    roiKeys = keys(roi)
)

    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@thread for i in eachindex(data)
        if passes_roi( data[i], roi, roiKeys)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end

    return ROIEfficiencyND(
        roi,  
        count[] / nTotalSim  # Extract final atomic count
    )
end

function get_roi_effciencyND(
    process::DataProcessND,
    roi::NamedTuple
)
    return get_effciencyND(
        process.data,
        roi,
        process.nTotalSim,
        process.varNames
    )
end

@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::NamedTuple,
    varNames
)
    @inbounds for n in varNames
        range = roi[n] 
        if !(range[1] ≤ data[n] < range[2])
            return false
        end
    end
    return true
end

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
        amount,
        keys(bins) # varNames
    )
end

function get_roi_bkg_counts(
    processes::Vector{DataProcessND}, 
    roi::NamedTuple,
)
    bkg = 0.0

    for p in processes
        p.signal && continue
        eff = get_roi_effciencyND(p, roi).eff
        bkg += p.amount * eff * p.activity * p.timeMeas
    end
    return bkg
end

function get_s_to_b(SNparams, 
    α, 
    processes::Vector{DataProcessND}, 
    ROIs::Vector{<:Real};
    approximate="table"
)
    # Construct ROI 
    roi = NamedTuple(
            k => (round(ROIs[i]), round(ROIs[i+1])) 
            for (i, k) in zip(1:2:length(processes[1].varNames)*2-1, processes[1].varNames)
        ) # creates a namedTuple in format (varName => (start, end), ...)

    signal_id = findfirst(p -> p.signal, processes)
    if signal_id === nothing
        @error "No signal process found!"
        return 0.0
    end
                    
    signal_process = processes[signal_id]
    ε = get_roi_effciencyND(signal_process, roi).eff
    ε == 0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return ε / S_b
end

function get_sensitivityND(
    SNparams::Dict, 
    α::Real, 
    processes::Vector{DataProcessND}, 
    roi::NamedTuple;
    approximate="table"
)

    #check that rois are within the range of the data
    for k in processes[1].varNames
        if (roi[k][1] < processes[1].bins[k][1] || roi[k][2] > processes[1].bins[k][2])
            @error "ROI out of range for $k"
        end
    end

    signal_id = findall([p.signal for p in processes])
    length(signal_id) > 1 && @error "Only one signal process allowed! Provided $(length(signal_id))"

    signal_process = processes[first(signal_id)]
    ε = get_roi_effciencyND(signal_process, roi).eff
    ε == 0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
    S_b = get_FC(b, α; approximate=approximate)

    tHalf = constantTerm * ε / S_b

    # return tHalf
    return SensitivityEstimateND(roi, tHalf, ε, b)
end


function set_activity!(process::DataProcessND, activity::Real)
    process.activity = activity
    return process
end

function set_timeMeas!(process::DataProcessND, timeMeas::Real)
    process.timeMeas = timeMeas
    return process
end

function set_nTotalSim!(process::DataProcessND, nTotalSim::Real)
    process.nTotalSim = nTotalSim
    return process
end


function set_amount!(process::DataProcessND, amount::Real)
    process.amount = amount
    return process
end
