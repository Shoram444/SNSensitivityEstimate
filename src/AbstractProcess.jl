abstract type AbstractProcess end

function get_process(isotopeName::String, processesVector::Vector{<:AbstractProcess})
    idx = findall( x -> x.isotopeName == isotopeName , processesVector )
    if( idx == nothing )
        @error "process $isotopeName not found!"
    end
    if( length(idx) > 1 )
        @info "process $isotopeName found $(length(idx)) times!"
    end
    
    # return all processes that were found with the same names
    return filter(x -> x.isotopeName == isotopeName, processesVector)
end

function set_signal!(process::AbstractProcess, signal::Bool)
    process.signal = signal
    return process
end

"""
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts(processes::AbstractProcess...)
    h2d = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    #check all processes have the same binning
    for p in processes
        if (p.bins != processes[1].bins)
            @error("get_bkg_counts(): passed processes have different binning!! Make sure they are the same.")
        end
    end

    for p in processes
        if (p.signal)
            @warn("get_bkg_counts(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            merge!(h2d, get_bkg_counts(p))
        end
    end
    return h2d
end


"""
!!
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts_1D(processes::Vector{<:AbstractProcess})
    h1d = Hist1D(; binedges=(processes[1].bins))

    for p in processes
        if (p.signal)
            @warn("get_bkg_counts(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            merge!(h1d, get_bkg_counts_1D(p))
        end
    end
    return h1d
end


"""
    Returns the number of expected counts of the given process. 
    (n = ε⋅t⋅m⋅A)
"""
function get_sig_counts(processes::Vector{<:AbstractProcess})
    h2d = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes
        if (!p.signal)
            @warn("get_sig_counts(): passed isotope $(p.isotopeName) is a background process!!")
        else
            merge!(h2d, get_sig_counts(p))
        end
    end
    return h2d
end


"""
    returns epsilon/(S(b)) 

"""
function get_epsilon_to_b(α, processes::Vector{<:AbstractProcess}; approximate="formula")
    ε = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes # sum efficiencies of the signal processes (does this make sense?)
        if(p.signal)
            ε += p.efficiency
        end
    end

    backgroundCounts = get_bkg_counts(processes)
    backgroundCounts.bincounts .= get_FC.( backgroundCounts.bincounts , α; approximate=approximate)
    epsToB = ε/(backgroundCounts)

    replace!(epsToB.bincounts, NaN => 0.0)

    return epsToB
end


function get_sToBRatio(processes::Vector{<:AbstractProcess})
    backgroundCounts = get_bkg_counts(processes)
    signalCounts = get_sig_counts(processes)

    if size(bincounts(backgroundCounts)) != size(bincounts(signalCounts))
        throw(ArgumentError("Input matrices must have the same dimensions"))
    end

    StoB = signalCounts / backgroundCounts
    replace!(StoB.bincounts, NaN => 0.0)

    return StoB
end




function get_tHalf_map(SNparams, α, processes::AbstractProcess...; approximate="formula")
    ε = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes # sum efficiencies of the signal processes (does this make sense?)
        if(p.signal)
            merge!(ε, p.efficiency)
        end
    end
    b = get_bkg_counts(processes...)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
     
    b.bincounts .= get_FC.(b.bincounts, α; approximate=approximate)

    tHalf = constantTerm * ε / b
    replace!(tHalf.bincounts, NaN => 0.0)

    return tHalf
end




"""
    get_estimated_bkg_counts(best_ROI, SNparams, processes::Vector{<:AbstractProcess})
==========
Outputs the estimated bkg counts for the given ROI. 
This process is obtained from looking up the `bkg_rate` at `best_ROI` and multiplying by mass and time.  
"""
function get_bkg_counts_ROI(minBinCenter, maxBinCenter,processes::AbstractProcess...)
    b = get_bkg_counts(processes...)

    return lookup(b , minBinCenter, maxBinCenter)
end

function get_bkg_counts_ROI(best_ROI, processes::AbstractProcess...)
    binStepHalf = step(processes[1].bins) / 2              # get the binning step and divide by half
    minBinCenter = best_ROI[:minBinEdge] + binStepHalf       # get the center of the minimal bin in ROI
    maxBinCenter = best_ROI[:maxBinEdge] - binStepHalf       # get the center of the maximal bin in ROI

    return get_bkg_counts_ROI(minBinCenter, maxBinCenter, processes...)
end

function FHist.lookup(process::AbstractProcess, best_ROI::Dict) 
    binStepHalf = step(process.bins) / 2              # get the binning step and divide by half
    minBinCenter = best_ROI[:minBinEdge] + binStepHalf       # get the center of the minimal bin in ROI
    maxBinCenter = best_ROI[:maxBinEdge] - binStepHalf       # get the center of the maximal bin in ROI
    if typeof(process) == HistProcess
        return lookup(process.counts, minBinCenter, maxBinCenter)
    end
    return lookup(process.efficiency, minBinCenter, maxBinCenter)
end

FHist.lookup(process::AbstractProcess, x::Real, y::Real) = typeof(process) == HistProcess ? lookup(process.counts, x, y) : lookup(process.efficiency, x, y)
DrWatson.default_allowed(::AbstractProcess) = (Real, String, Bool, AbstractRange)
DrWatson.allaccess(::AbstractProcess) = (:isotopeName, :signal, :bins, :activity, :nTotalSim, :amount)


