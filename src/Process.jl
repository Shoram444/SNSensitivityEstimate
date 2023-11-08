mutable struct Process
    dataVector
    isotopeName
    signal
    activity
    timeMeas
    nTotalSim
    bins
    efficiency
end

function Process(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins)

    eff = get_efficiency(dataVector, bins, nTotalSim)
    return Process(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins,eff)
end

function plot_efficiency(process::Process, bins; nTotalSim = 1e8, kwargs...)
    n = nTotalSim
    h1d = StatsBase.fit(Histogram{Float64}, process.dataVector, bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins,bins)) #prepare empty 2d Histogram
    h2d.weights = get_efficiency(process)
    plot(h2d; kwargs...)
end

function plot_efficiency(effMatrix::Matrix, bins, nTotalSim = 1e8; kwargs...)
    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins,bins)) #prepare empty 2d Histogram
    h2d.weights = get_efficiency(effMatrix, bins, nTotalSim)
    plot(h2d; kwargs...)
end


function get_efficiency(process::Process, bins, nTotalSim = 1e8; kwargs...)
    n = nTotalSim
    h1d = StatsBase.fit(Histogram{Float64}, process.dataVector, bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins,bins)) #prepare empty 2d Histogram

    minBin, stepBin, maxBin = bins[1], step(bins), bins[end-1]

    for minVal in minBin:stepBin:maxBin-stepBin
        for maxVal in minVal+stepBin:stepBin:maxBin
            minBinID = FHist._edge_binindex(bins, minVal)
            maxBinID = FHist._edge_binindex(bins, maxVal)
            binSum = sum(h1d.weights[Int(minBinID):Int(maxBinID)]) # integrate the histogram from minimal bin to maximal

            push!(h2d, (minVal, maxVal), binSum/n)
        end
    end
    h2d.weights
end

function get_efficiency(dataVector, bins, nTotalSim; kwargs...)
    n = nTotalSim
    h1d = StatsBase.fit(Histogram{Float64}, dataVector, bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins,bins)) #prepare empty 2d Histogram

    minBin, stepBin, maxBin = bins[1], step(bins), bins[end-1]

    for minVal in minBin:stepBin:maxBin-stepBin
        for maxVal in minVal+stepBin:stepBin:maxBin
            minBinID = FHist._edge_binindex(bins, minVal)
            maxBinID = FHist._edge_binindex(bins, maxVal)
            binSum = sum(h1d.weights[Int(minBinID):Int(maxBinID)]) # integrate the histogram from minimal bin to maximal

            push!(h2d, (minVal, maxVal), binSum/n)
        end
    end
    h2d.weights
end

function get_efficiency(process::Process; kwargs...)
    return process.efficiency 
end

function get_efficiency(effMatrix::Matrix, bins, nTotalSim = 1e8; kwargs...)
    n = nTotalSim

    return effMatrix./n
end

function get_nPassed(process::Process, bins; kwargs...)
    h1d = StatsBase.fit(Histogram{Float64}, process.dataVector, bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins,bins)) #prepare empty 2d Histogram

    minBin, stepBin, maxBin = bins[1], step(bins), bins[end-1]

    for minVal in minBin:stepBin:maxBin-stepBin
        for maxVal in minVal+stepBin:stepBin:maxBin
            minBinID = FHist._edge_binindex(bins, minVal)
            maxBinID = FHist._edge_binindex(bins, maxVal)
            binSum = sum(h1d.weights[Int(minBinID):Int(maxBinID)]) # integrate the histogram from minimal bin to maximal

            push!(h2d, (minVal, maxVal), binSum)
        end
    end
    h2d.weights
end

function get_nPassed(process::Process; kwargs...)
    h1d = StatsBase.fit(Histogram{Float64}, process.dataVector, process.bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (process.bins,process.bins)) #prepare empty 2d Histogram

    minBin, stepBin, maxBin = process.bins[1], step(process.bins), process.bins[end-1]

    for minVal in minBin:stepBin:maxBin-stepBin
        for maxVal in minVal+stepBin:stepBin:maxBin
            minBinID = FHist._edge_binindex(process.bins, minVal)
            maxBinID = FHist._edge_binindex(process.bins, maxVal)
            binSum = sum(h1d.weights[Int(minBinID):Int(maxBinID)]) # integrate the histogram from minimal bin to maximal

            push!(h2d, (minVal, maxVal), binSum)
        end
    end
    h2d.weights
end

function plot_nPassed(process::Process, bins; kwargs...)
    h1d = StatsBase.fit(Histogram{Float64}, process.dataVector, bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins,bins)) #prepare empty 2d Histogram

    minBin, stepBin, maxBin = bins[1], step(bins), bins[end-1]

    for minVal in minBin:stepBin:maxBin-stepBin
        for maxVal in minVal+stepBin:stepBin:maxBin
            minBinID = FHist._edge_binindex(bins, minVal)
            maxBinID = FHist._edge_binindex(bins, maxVal)
            binSum = sum(h1d.weights[Int(minBinID):Int(maxBinID)]) # integrate the histogram from minimal bin to maximal

            push!(h2d, (minVal, maxVal), binSum)
        end
    end
    plot(h2d; kwargs...)
end

function plot_nPassed(process::Process; kwargs...)
    h1d = StatsBase.fit(Histogram{Float64}, process.dataVector, process.bins)

    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (process.bins,process.bins)) #prepare empty 2d Histogram

    minBin, stepBin, maxBin = process.bins[1], step(process.bins), process.bins[end-1]

    for minVal in minBin:stepBin:maxBin-stepBin
        for maxVal in minVal+stepBin:stepBin:maxBin
            minBinID = FHist._edge_binindex(process.bins, minVal)
            maxBinID = FHist._edge_binindex(process.bins, maxVal)
            binSum = sum(h1d.weights[Int(minBinID):Int(maxBinID)]) # integrate the histogram from minimal bin to maximal

            push!(h2d, (minVal, maxVal), binSum)
        end
    end
    plot(h2d; kwargs...)
end

function lookup_efficiency(process::Process, x, y)
    effMatrix = process.efficiency
    binning = process.bins
    rx, ry = binning, binning
    !(first(rx) <= x <= last(rx)) && return missing
    !(first(ry) <= y <= last(ry)) && return missing
    @show binning[FHist._edge_binindex(rx, x)], binning[FHist._edge_binindex(ry, y)]
    @show FHist._edge_binindex(rx, x), FHist._edge_binindex(ry, y)
    return effMatrix[FHist._edge_binindex(rx, x), FHist._edge_binindex(ry, y)]
end

function lookup_efficiency(h::Histogram, x, y)
    effMatrix = h.weights
    binning = h.edges[1]
    rx, ry = binning, binning
    !(first(rx) <= x <= last(rx)) && return missing
    !(first(ry) <= y <= last(ry)) && return missing
    @show binning[FHist._edge_binindex(rx, x)], binning[FHist._edge_binindex(ry, y)]
    @show FHist._edge_binindex(rx, x), FHist._edge_binindex(ry, y)
    return effMatrix[FHist._edge_binindex(rx, x), FHist._edge_binindex(ry, y)]
end


function get_bkg_counts(processes::Process...)

    backgroundCounts = zeros(size(processes[1].efficiency)) # initiaite empty Matrix
    for p in processes
        if(p.signal)
            @warn("get_bkg_counts(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            @show p.efficiency[1,:]
            backgroundCounts .+= p.activity* p.timeMeas *p.efficiency
            @show backgroundCounts[1,:]
        end
    end
    return backgroundCounts
end

function get_sig_counts(processes::Process...)

    signalCounts = zeros(size(processes[1].efficiency)) # initiaite empty Matrix
    for p in processes
        if(!p.signal)
            @warn("get_sig_counts(): passed isotope $(p.isotopeName) is a background process!!")
        else
            signalCounts .+= p.activity* p.efficiency
        end
    end
    return signalCounts
end

function get_sToBRatio(processes::Process...)
    backgroundCounts = get_bkg_counts(processes...)
    signalCounts = get_sig_counts(processes...)

    if size(backgroundCounts) != size(signalCounts)
        throw(ArgumentError("Input matrices must have the same dimensions"))
    end

    sToB = zeros(size(signalCounts)) # initialize empty
    for i in eachindex(signalCounts)
        !(backgroundCounts[i]==0.0) && (sToB[i] = signalCounts[i] / backgroundCounts[i]) # in each cell divide signal by bkg if bkg ≠ 0
    end
    
    return sToB
end

function plot_sToBRatio(processes::Process...; kwargs...)
    bins = processes[1].bins
    h2d = StatsBase.fit(Histogram{Float64}, (Float64[], Float64[]), (bins, bins))
    
    h2d.weights = get_sToBRatio(processes...)
    return plot(h2d; kwargs...)
end