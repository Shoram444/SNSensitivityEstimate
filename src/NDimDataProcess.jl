import Base: *, maximum, println

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

Base.println(s::SensitivityEstimateND) = println(
    "Sensitivity estimate:
     ROI:   $(s.roi)
     T12 ≥ $(s.tHalf) yr
     ε = $(s.signalEff)
     b = $(s.bkgCounts)
    "
)

"""
    Holds information about ROI efficiency
"""
mutable struct ROIEfficiencyND
    roi::Vector{<:Real}
    eff
end

*(x::Float64, e::ROIEfficiencyND) = x * e.eff

mutable struct DataProcessND <: AbstractProcess 
    data::LazyTree
    isotopeName::String
    signal::Bool
    activity::Union{Real, Measurement}
    timeMeas::Real
    nTotalSim::Real
    bins::NamedTuple
    amount::Real
    varNames::Vector{String}
    varIdxs::Vector{Int}
end



function get_roi_effciencyND(
    data::LazyTree, 
    roi::NamedTuple, 
    nTotalSim::Real,
    roiKeys = keys(roi)
)

    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        if passes_roi( data[i], roi, roiKeys)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end

    ROI = Float64[]
    for k in keys(roi)
        push!(ROI, roi[k][1])
        push!(ROI, roi[k][2])
    end

    return ROIEfficiencyND(
        ROI,  
        count[] / nTotalSim  # Extract final atomic count
    )
end

function get_roi_effciencyND(
    data::LazyTree, 
    roi::NamedTuple, 
    nTotalSim::Real,
    roiKeys::Vector{String}
)

    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        if passes_roi( data[i], roi, roiKeys)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end

    ROI = Float64[]
    for k in keys(roi)
        push!(ROI, roi[k][1])
        push!(ROI, roi[k][2])
    end

    return ROIEfficiencyND(
        ROI,  
        count[] / nTotalSim  # Extract final atomic count
    )
end

function get_roi_effciencyND(
    data::LazyTree, 
    roi::Vector{<:Real}, 
    nTotalSim::Real,
    varIdxs::Vector{Int}
)

    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        if passes_roi( data[i], roi, varIdxs)
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
    return get_roi_effciencyND(
        process.data,
        roi,
        process.nTotalSim,
        process.varNames
    )
end

function get_roi_effciencyND(
    process::DataProcessND,
    roi::Vector{<:Real},
)
    return get_roi_effciencyND(
        process.data,
        roi,
        process.nTotalSim,
        process.varIdxs
    )
end

@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::Vector{<:Real},
    varIdxs::Vector{Int}
)
    @inbounds for (r,i) in zip(1:2:length(roi)*2-1, varIdxs)
        if !(roi[r] ≤ data[i] < roi[r+1])
            return false
        end
    end
    return true
end



@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::Vector{<:Real},
    varNames::Vector{Symbol}
)
    @inbounds for (i,d) in zip(1:2:length(roi), data)
        if !(roi[i] ≤ d < roi[i+1])
            return false
        end
    end
    return true
end


@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::NamedTuple,
    varNames::Vector{String}
)
    @inbounds for n in varNames
        range = roi[Symbol(n)]
        d = data[Symbol(n)] 
        if !(range[1] ≤ d < range[2])
            return false
        end
    end
    return true
end

function DataProcessND(
    data::LazyTree,
    binsTuple::NamedTuple,
    varNames::Vector{String},
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
        binsTuple,
        amount,
        varNames,
        [findfirst(x -> x == n, names(data)) for n in varNames]
    )
end

function get_roi_bkg_counts(
    processes::Vector{DataProcessND}, 
    roi::NamedTuple,
)
    bkg = zero(Float64)
    ε = zero(Float64)

    for p in processes
        p.signal && continue
        ε = get_roi_effciencyND(p, roi).eff
        bkg += p.amount * ε * p.activity * p.timeMeas
    end
    return bkg
end

function get_roi_bkg_counts(
    processes::Vector{DataProcessND}, 
    roi::Vector{<:Real},
)
    bkg = zero(Float64)
    ε = zero(Float64)

    for p in processes
        p.signal && continue
        ε = get_roi_effciencyND(p, roi).eff
        bkg += p.amount * ε * p.activity * p.timeMeas
    end
    return bkg
end

function get_s_to_b(
    SNparams::Dict, 
    α::Float64, 
    processes::Vector{DataProcessND}, 
    roi::Vector{<:Real};
    approximate="table"
)
    signal_id = findfirst(p -> p.signal, processes)
    if signal_id === nothing
        @error "No signal process found!"
        return 0.0
    end
                    
    ε = get_roi_effciencyND(processes[signal_id], roi).eff
    ε == 0.0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_roi_bkg_counts(processes, roi)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return ε / S_b
end

# function get_s_to_b(
#     SNparams::Dict, 
#     α::Float64, 
#     processes::Vector{DataProcessND}, 
#     ROIs::NamedTuple;
#     approximate="table"
# )
#     # Construct ROI 
#     roi = NamedTuple(
#             Symbol(k) => (round(ROIs[k][i]), round(ROIs[k][i+1])) 
#             for (i, k) in zip(1:2:length(processes[1].varNames)*2-1, processes[1].varNames)
#         ) # creates a namedTuple in format (varName => (start, end), ...)

#     signal_id = findfirst(p -> p.signal, processes)
#     if signal_id === nothing
#         @error "No signal process found!"
#         return 0.0
#     end
                    
#     # signal_process = processes[signal_id]
#     @time ε = get_roi_effciencyND(processes[signal_id], roi).eff
#     ε == 0.0 && return 0.0 # If efficiency is zero, return zero sensitivity

#     @time b = get_roi_bkg_counts(processes, roi)
#     @unpack W, foilMass, Nₐ, tYear, a = SNparams
#     S_b = get_FC(b, α; approximate=approximate)

#     return ε / S_b
# end

function get_sensitivityND(
    SNparams::Dict, 
    α::Real, 
    processes::Vector{DataProcessND}, 
    roi::NamedTuple;
    approximate="table"
)

    #check that rois are within the range of the data
    for k in processes[1].varNames
        if (roi[Symbol(k)][1] < processes[1].bins[Symbol(k)][1] || roi[Symbol(k)][2] > processes[1].bins[Symbol(k)][2])
            @error "ROI out of range for $k"
        end
    end

    signal_id = findall([p.signal for p in processes])
    if length(signal_id) > 1 
        @error "Only one signal process allowed! Provided $(length(signal_id))" 
        return SensitivityEstimateND(roi, 0.0, 0.0, 0.0)
    end

    signal_process = processes[first(signal_id)]
    ε = get_roi_effciencyND(signal_process, roi).eff
    ε == 0 && return SensitivityEstimateND(roi, 0.0, 0.0, 0.0) # If efficiency is zero, return zero sensitivity

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
