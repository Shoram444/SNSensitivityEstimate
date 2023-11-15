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


function get_efficiency(dataVector, bins, nTotalSim; kwargs...)
    n = nTotalSim
    
    h2d = get_nPassed(dataVector, bins) * (1/n) # Matrix - scalar multiplication
end


function get_nPassed(process::Process, bins)
    return get_nPassed(process.dataVector, bins)
end

function get_nPassed(dataVector, bins)
    h1d = Hist1D(dataVector, bins )
    h2d = Hist2D(Float64, bins = (bins, bins)  )#prepare empty 2d Histogram

    binCenters = collect(bincenters(h1d))
    for (i, xVal) in enumerate(binCenters)
        for (j, yVal) in enumerate(binCenters[i+1:end])
            j += i

            integral = sum(lookup.(h1d, binCenters[i:j]))#xVal:yVal))  
            push!(h2d, xVal, yVal ,integral)
        end
    end

    h2d
end


FHist.lookup(process::Process, x, y) = lookup(process.efficiency, x, y)



function get_bkg_rate(processes::Process...)
    h2d = Hist2D(Float64, bins = (processes[1].bins, processes[1].bins)  )

    for p in processes
        if(p.signal)
            @warn("get_bkg_rate(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            h2d += p.efficiency * (p.activity* p.timeMeas )
        end
    end
    return h2d
end

function get_sig_rate(processes::Process...)
    h2d = Hist2D(Float64, bins = (processes[1].bins, processes[1].bins)  )

    for p in processes
        if(!p.signal)
            @warn("get_sig_rate(): passed isotope $(p.isotopeName) is a background process!!")
        else
            h2d += p.efficiency * p.activity 
        end
    end
    return h2d
end

function get_sToBRatio(processes::Process...)
    backgroundCounts = get_bkg_rate(processes...)
    signalCounts = get_sig_rate(processes...)

    if size(bincounts(backgroundCounts)) != size(bincounts(signalCounts))
        throw(ArgumentError("Input matrices must have the same dimensions"))
    end
    
    StoB = signalCounts / backgroundCounts
    replace!(StoB.hist.weights, NaN => 0.0)

    return StoB
end

