import Base: *, maximum, println

"""
    Holds information about sensitivity estimate
    Fields:
    - `roi`: NamedTuple of the ROI
    - `tHalf`: sesnitivity estimate
    - `signalEff`: Signal efficiency
    - `bkgCounts`: Background counts
"""
struct SensitivityEstimateND{T, R, M}
    roi::R
    tHalf::T
    signalEff::T
    bkgCounts::M
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
mutable struct ROIEfficiencyND{T, M}
    roi::Vector{T}
    eff::M
end

*(x::Float64, e::ROIEfficiencyND) = x * e.eff

mutable struct DataProcessND{T, B, NT, S, I} <: AbstractProcess 
    # data::LazyTree
    data::Vector{UnROOT.LazyEvent}
    isotopeName::S
    signal::B
    activity::Union{T, Measurement}
    timeMeas::T
    nTotalSim::T
    bins::NT
    amount::T
    varNames::Vector{S}
    varIdxs::Vector{I}
end



function get_roi_effciencyND(
    # data::LazyTree, 
    data::Vector{UnROOT.LazyEvent}, 
    roi::NamedTuple, 
    nTotalSim::Real,
    roiKeys = keys(roi)
)
    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        evt = data[i]
        if passes_roi( evt, roi, roiKeys)
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
    data::Vector{UnROOT.LazyEvent}, 
    # data::LazyTree, 
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
    data::Vector{UnROOT.LazyEvent}, 
    # data::LazyTree, 
    roi::Vector{<:Real}, 
    nTotalSim::Real,
    varIdxs::Vector{Int}
)
    count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

    Threads.@threads for i in eachindex(data)
        if passes_roi(data[i], roi, varIdxs)
            Threads.atomic_add!(count, 1)  # Thread-safe increment
        end
    end
    eff = count[] / nTotalSim  # Extract final atomic count
    return ROIEfficiencyND(
        roi,  
        eff  # Extract final atomic count
    )
end

# function get_roi_effciencyND(
#     data::Vector{UnROOT.LazyEvent}, 
#     roi::Vector{<:Real}, 
#     nTotalSim::Real,
#     varIdxs::Vector{Int}
# )

#     count = Threads.Atomic{Int}(0)  # Atomic counter for thread-safe increment

#     Threads.@threads for i in eachindex(data)
#         if passes_roi( data[i], roi, varIdxs)
#             Threads.atomic_add!(count, 1)  # Thread-safe increment
#         end
#     end

#     return ROIEfficiencyND(
#         roi,  
#         count[] / nTotalSim  # Extract final atomic count
#     )
# end

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
    roi::Vector{<:Real}
)
    return get_roi_effciencyND(
        getproperty(process, :data),
        roi,
        getproperty(process, :nTotalSim),
        getproperty(process, :varIdxs)
    )
end

@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::Vector{<:Real},
    varIdxs::Vector{Int}
)
    @inbounds for (r, i) in zip(1:2:length(roi)*2-1, varIdxs)
        if !(roi[r] ≤ abs(data[i]) < roi[r+1])
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
    i = 0
    @inbounds for (i,d) in zip(1:2:length(roi), data)
        if !(roi[i] ≤ abs(d) < roi[i+1])
            return false
        end
        i += 1
    end
    i == length(varNames) && return false
    return true
end


@inline function passes_roi(
    data::UnROOT.LazyEvent,  
    roi::NamedTuple,
    varNames::Vector{String}
)
    i = 1
    @inbounds for n in varNames
        range = roi[Symbol(n)]
        d = data[Symbol(n)] 
        if !(range[1] ≤ abs(d) < range[2])
            return false
        end
        i += 1
    end
    i == length(varNames) && return false
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

    collected_data = collect(data)
    varIdxs = [findfirst(x -> x == n, names(data)) for n in varNames]

    T = promote_type(typeof(activity), typeof(timeMeas), typeof(nTotalSim), typeof(amount))
    I = promote_type(Int64, typeof(varIdxs[1]))

    return DataProcessND{T, Bool, typeof(binsTuple), String, I}(
        collected_data,
        isotopeName,
        signal,
        activity,
        timeMeas,
        nTotalSim,
        binsTuple,
        amount,
        varNames,
        varIdxs
    )
end


function get_roi_bkg_counts(
    processes::Vector{<:DataProcessND}, 
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
    processes::Vector{<:DataProcessND}, 
    roi::Vector{<:Real},
)
    bkg = zero(Float64)

    for p in processes
        p.signal && continue
        bkg += p.amount * get_roi_effciencyND(p, roi).eff * p.activity * p.timeMeas
    end
    return bkg
end

function get_roi_signal_efficiency(
    processes::Vector{<:DataProcessND}, 
    roi::Vector{<:Real},
)
    ε = zero(Float64)

    for p in processes
        p.signal == false && continue
        ε += get_roi_effciencyND(p, roi).eff
    end
    return ε
end

function get_s_to_b(
    SNparams::Dict, 
    α::Float64, 
    processes::Vector{<:DataProcessND}, 
    roi::Vector{Float64};
    approximate="table"
)
    for i in 1:2:length(roi)-1
        if roi[i] >= roi[i+1]
            return -1e20  # penalty for invalid ROI
        end
    end

    signal_id = findfirst(p -> p.signal, processes)
    if signal_id === nothing
        @error "No signal process found!"
        return -1e20 # penalty for no signal process
    end
                    
    ε = get_roi_signal_efficiency(processes, roi)
    ε == 0.0 && return -1e20 # If efficiency is zero, penalize optimization  

    b = get_roi_bkg_counts(processes, roi)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return Measurements.value(ε / S_b)
end

function get_sensitivityND(
    SNparams::Dict, 
    α::Real, 
    processes::Vector{<:DataProcessND}, 
    roi::NamedTuple;
    approximate="table",
    add_mock_bkg =0.0
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

    b = get_roi_bkg_counts(processes, roi) + add_mock_bkg

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
