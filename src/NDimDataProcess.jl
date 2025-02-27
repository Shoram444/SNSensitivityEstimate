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
    nthreads = Threads.nthreads()  # Get number of available threads
    local_counts = fill(0, nthreads)  # Create thread-local counters

    Threads.@threads for tid in 1:nthreads
        local_count = 0  # Local variable for counting
        @inbounds for i in tid:nthreads:length(dataAngle)  # Chunking work by thread ID
            if (
                binsESum[1] ≤ dataESum[i] < binsESum[2] &&
                binsESingle[1] ≤ dataESingle[i] < binsESingle[2] &&
                binsAngle[1] ≤ dataAngle[i] < binsAngle[2]
            )
                local_count += 1  # Update local count
            end
        end
        local_counts[tid] = local_count  # Store count from this thread
    end

    return sum(local_counts) / nTotalSim  # Sum up all thread-local counts
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
    println("creating process: $isotopeName")
    # Generate all possible region of interest (ROI) combinations
    roi_combinations = collect(Iterators.product(
        make_ROI_combinations(binsAngle),
        make_ROI_combinations(binsESingle),
        make_ROI_combinations(binsESum)
    ))

    # Ensure data vectors have the same length
    (length(dataAngle) == length(dataESingle) == length(dataESum)) || 
        @error "Data vectors must have the same length!"

    # Initialize multi-threading
    counts = zeros(Int, length(roi_combinations))

    # Threaded loop for counting events in each ROI
    @inbounds for event_idx in 1:length(dataAngle)
        event_idx % 1_000_000 == 0 && println("$event_idx/$(length(dataAngle)) events processed!")

        angle_val = dataAngle[event_idx]
        esingle_val = dataESingle[event_idx]
        esum_val = dataESum[event_idx]

        @inbounds for (roi_idx, (angle_roi, esingle_roi, esum_roi)) in enumerate(roi_combinations)

            if (passes_roi(angle_val, esingle_val, esum_val, angle_roi, esingle_roi, esum_roi))
                counts[roi_idx] += 1
            end
        end
    end

    # Compute efficiency for each ROI
    efficiencies = [NDimROIEfficiency(angle, esingle, esum, count / nTotalSim) 
                    for ((angle, esingle, esum), count) in zip(roi_combinations, counts)]

    return NDimDataProcess(dataAngle, dataESingle, dataESum, isotopeName, signal, activity, timeMeas, nTotalSim, binsAngle, binsESingle, binsESum, amount, efficiencies)
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

function precompute_histogram(dataAngle, dataESingle, dataESum, binsAngle, binsESingle, binsESum)
    hist = fit(Histogram, (dataAngle, dataESingle, dataESum), (binsAngle, binsESingle, binsESum))
    return hist.weights  # 3D array storing counts per bin
end

function get_nDim_effciency(histogram::Array, bin_idx_angle, bin_idx_esingle, bin_idx_esum, nTotalSim)
    return histogram[bin_idx_angle, bin_idx_esingle, bin_idx_esum] / nTotalSim
end
