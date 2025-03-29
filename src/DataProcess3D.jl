import Base: *, maximum

"""
    aa
"""
mutable struct SensitivityEstimate3D
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
mutable struct ROIEfficiency3D
    binESingle
    binESum
    binPhi
    eff
end


*(x::Float64, e::ROIEfficiency3D) = x * e.eff

mutable struct DataProcess3D <: AbstractProcess
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
end



function get_effciency3D(
        dataAngle::Vector{<:Real}, 
        dataESingle::Vector{<:Real}, 
        dataESum::Vector{<:Real}, 
        binsAngle::Vector{<:Real}, 
        binsESingle::Vector{<:Real}, 
        binsESum::Vector{<:Real}, 
        nTotalSim::Real
    )
    nthreads = Threads.nthreads()  # Get number of available threads
    local_counts = fill(0, nthreads)  # Create thread-local counters

    Threads.@threads for tid in 1:nthreads
        local_count = 0  # Local variable for counting
        @inbounds for i in tid:nthreads:length(dataAngle)  # Chunking work by thread ID
            if (
                passes_roi(
                    dataAngle[i], 
                    dataESingle[i], 
                    dataESum[i], 
                    binsAngle, 
                    binsESingle, 
                    binsESum
                )
            )
                local_count += 1  # Update local count
            end
        end
        local_counts[tid] = local_count  # Store count from this thread
    end
    
    return ROIEfficiency3D(
        binsESingle,
        binsESum,
        binsAngle,
        sum(local_counts) / nTotalSim
    )
end

function get_effciency3D(
    process::DataProcess3D,
    binsAngle::Vector{<:Real}, 
    binsESingle::Vector{<:Real}, 
    binsESum::Vector{<:Real}, 
)
    return get_effciency3D(
        process.dataAngle,
        process.dataESingle,
        process.dataESum,
        binsAngle,
        binsESingle,
        binsESum,
        process.nTotalSim
    )
end



@inline function passes_roi(
        angle_val::Real, 
        esingle_val::Real, 
        esum_val::Real, 
        angle_roi::Vector{<:Real}, 
        esingle_roi::Vector{<:Real}, 
        esum_roi::Vector{<:Real}
    )
    return (
        esum_roi[1] ≤ esum_val < esum_roi[2] &&
        angle_roi[1] ≤ angle_val < angle_roi[2] &&
        esingle_roi[1] ≤ esingle_val < esingle_roi[2]
        
    )
end

function DataProcess3D(
        dataAngle::Vector{<:Real},
        dataESingle::Vector{<:Real},
        dataESum::Vector{<:Real},
        binsAngle,
        binsESingle,
        binsESum,
        processDict::Dict
    )
    @unpack isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount = processDict
    println("creating process: $isotopeName")

    return DataProcess3D(
        dataAngle,
        dataESingle,
        dataESum,
        isotopeName,
        signal,
        activity,
        timeMeas,
        nTotalSim,
        binsAngle,
        binsESingle,
        binsESum,
        amount
    )
end

function get_bkg_counts(
        processes::Vector{DataProcess3D}, 
        binsAngle::Vector{<:Real}, 
        binsESingle::Vector{<:Real}, 
        binsESum::Vector{<:Real}
    )
    bkg = 0.0

    for p in processes
        p.signal && continue
        eff = get_effciency3D(p.dataAngle, p.dataESingle, p.dataESum, binsAngle, binsESingle, binsESum, p.nTotalSim).eff
        bkg += p.amount * eff * p.activity * p.timeMeas
    end
    return bkg
end

function get_s_to_b(SNparams, 
    α, 
    processes::Vector{DataProcess3D}, 
    ROIs::Vector{<:Real};
    approximate="table"
)
    binsAngle = round.(ROIs[1:2], digits = -1)
    binsESingle = round.(ROIs[3:4], digits = -2)
    binsESum = round.(ROIs[5:6], digits = -2)

    #check that bins are within the range of the data
    if( binsAngle[1] < processes[1].binsAngle[1] || binsAngle[2] > processes[1].binsAngle[end] )
        return -1e6
    end
    if( binsESingle[1] < processes[1].binsESingle[1] || binsESingle[2] > processes[1].binsESingle[end] )
        return -1e6
    end
    if( binsESum[1] < processes[1].binsESum[1] || binsESum[2] > processes[1].binsESum[end] )
        return -1e6
    end


    signal_id = findall([p.signal for p in processes])
    length(signal_id) > 1 && @error "Only one signal process allowed! Provided $(length(signal_id))"
    
    signal_process = processes[first(signal_id)]
    ε = get_effciency3D(signal_process.dataAngle, signal_process.dataESingle, signal_process.dataESum, binsAngle, binsESingle, binsESum, signal_process.nTotalSim).eff
    ε == 0 && return -1e6 # If efficiency is zero, penalize with a large number

    b = get_bkg_counts(processes, binsAngle, binsESingle, binsESum)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    S_b = get_FC(b, α; approximate=approximate)

    return ε / S_b
end

function get_sensitivity3D(
        SNparams, 
        α, 
        processes::Vector{DataProcess3D}, 
        ROIs::Vector{<:Real};
        approximate="table"
    )

    binsAngle = ROIs[1:2]
    binsESingle = ROIs[3:4]
    binsESum = ROIs[5:6]

    #check that bins are within the range of the data
    if( binsAngle[1] < processes[1].binsAngle[1] || binsAngle[2] > processes[1].binsAngle[end] )
        return 0.0
    end
    if( binsESingle[1] < processes[1].binsESingle[1] || binsESingle[2] > processes[1].binsESingle[end] )
        return 0.0
    end
    if( binsESum[1] < processes[1].binsESum[1] || binsESum[2] > processes[1].binsESum[end] )
        return 0.0
    end


    signal_id = findall([p.signal for p in processes])
    length(signal_id) > 1 && @error "Only one signal process allowed! Provided $(length(signal_id))"
    
    signal_process = processes[first(signal_id)]
    ε = get_effciency3D(signal_process.dataAngle, signal_process.dataESingle, signal_process.dataESum, binsAngle, binsESingle, binsESum, signal_process.nTotalSim).eff
    ε == 0 && return 0.0 # If efficiency is zero, return zero sensitivity

    b = get_bkg_counts(processes, binsAngle, binsESingle, binsESum)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
    S_b = get_FC(b, α; approximate=approximate)

    tHalf = constantTerm * ε / S_b

    # return tHalf
    return SensitivityEstimate3D(
            binsAngle, 
            binsESingle, 
            binsESum, 
            tHalf,
            ε,
            b
        )
    
end

function get_sensitivity3D(
    SNparams, 
    α, 
    processes::Vector{DataProcess3D}, 
    ROIs::Dict;
    approximate="table"
)

    binsAngle = ROIs[:angle]
    binsESingle = ROIs[:esingle]
    binsESum = ROIs[:esum]

    return get_sensitivity3D(
        SNparams, 
        α, 
        processes, 
        [binsAngle[1], binsAngle[2], binsESingle[1], binsESingle[2], binsESum[1], binsESum[2]];
        approximate=approximate
    )
end


function set_activity!(process::DataProcess3D, activity::Real)
    process.activity = activity
    return process
end

function set_timeMeas!(process::DataProcess3D, timeMeas::Real)
    process.timeMeas = timeMeas
    return process
end

function set_nTotalSim!(process::DataProcess3D, nTotalSim::Real)
    process.nTotalSim = nTotalSim
    return process
end



function set_amount!(process::DataProcess3D, amount::Real)
    process.amount = amount
    return process
end
