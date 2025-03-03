"""
    struct HistProcess -> `(hist1D, histName, isotopeName, signal, activity, timeMeas, nTotalSim, bins, efficiency)`

The following fields are defined:
    + hist::Hist1D - histogram of the data obtained from root file
    + isotopeName::String - isotope name 
    + histName::String - name of the histogram in the root file
    + signal::Bool - signal or background
    + counts::Hist2D - the 2D histogram with bkgCounts per bin 
"""
mutable struct HistProcess <: AbstractProcess
    hist::Hist1D
    isotopeName::String
    histName::String
    signal::Bool
    bins::AbstractRange
    counts::Hist2D
end

function HistProcess(
    hist::Hist1D,
    isotopeName::String,
    histName::String,
    signal::Bool,
    bins::AbstractRange
)
    counts = get_counts(hist) 
    return HistProcess(hist, isotopeName, histName, signal, bins, counts)
end

function HistProcess(
    hist::Hist1D,
    histName::String,
    processDict::Dict
)
    @unpack isotopeName, signal, bins = processDict


    # a hack to change from MeV to keV
    if( binedges(hist)[end] < 100 )
        hist = Hist1D(; binedges = binedges(hist) .* 1000, counttype = Float64, sumw2 = hist.sumw2, bincounts = bincounts(hist))
    end

    hist = restrict(hist, bins[1], bins[end])
    if (binedges(hist) != bins)
        @error("HistProcess(): binning from Params.jl doesn't match binning of the passed histogram!! Make sure they are the same. 
                Params binning: $(bins), passed histogram binning: $(binedges(hist))")
    end
    counts = get_counts(hist) 
    return HistProcess(hist, isotopeName, histName, signal, bins, counts)
end



"""
    get_counts(hist::Hist1D)

    returns a 2D histgoram of type Hist2D with counts per bin.

"""
function get_counts(hist::Hist1D)
    h2d=get_nPassed(hist) # nPassed scaled by total simulated events
    return h2d
end

"""
!!
    get_nPassed(process::DataProcess, bins)

    returns a 2D histgoram of type Hist2D with passed events per bin.

"""
function get_nPassed(process::HistProcess)
    return get_nPassed(process.hist)
end


function get_nPassed(hist::Hist1D)
    bins = binedges(hist)
    h2d = Hist2D(; counttype=Float64, binedges=(bins, bins))#prepare empty 2d Histogram

    bc=bincenters(hist)
    for (i,b) in enumerate(bc)
        cs=cumsum(bincounts(hist)[i:end]) 
        push!.(h2d, b, bc[i:end], cs)
    end

    h2d
end


"""
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts(process::HistProcess)
    return process.counts
end

"""
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts_1D(process::HistProcess)
    return process.hist
end



"""
    get_isotope_details( process::HistProcess )

    returns a tuple of isotope details: (a, m, t, ε, nExpTot, nTotSim)

    Useful for printing and estimating expected counts.
"""
function get_isotope_details( process::HistProcess )
    a = 0
    m = 0
    t = 0
    nTotSim = 0
    ε = 0

    nExpTot = sum(bincounts(process.hist))

    return a, m, t, ε, nExpTot, nTotSim 
end

"""
    print_isotope_details( processes::HistProcess... )

    returns a string with basic isotope details.
"""
function print_isotope_details( process::HistProcess... )
    (a, m, t, ε, nExpTot, nTotSim) = get_isotope_details(process...) 

    a = round(a, sigdigits=3)
    ε = round(ε*100, sigdigits=3)
    m = round(m, sigdigits=3)
    nExpTot = round(nExpTot, sigdigits=3)
    isotope = process.isotopeName
    
    "|$isotope | $nExpTot | $ε | $a | $m |"
end

function get_process(isotopeName::String, processesVector::Vector{HistProcess}, histName::String)
    idx = findfirst( x -> (x.isotopeName == isotopeName && x.histName == histName), processesVector )
    if( idx == nothing )
        @error "process $isotopeName not found!"
    end
    return processesVector[idx]
end

