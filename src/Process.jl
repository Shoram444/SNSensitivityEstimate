"""
    mutable struct Process -> `(dataVector, isotopeName, signa, activit, timeMea, nTotalSi, bin, efficiency)`

The following fields are defined:
    + dataVector::Vector{<:Real} - vector of initial data points
    + isotopeName::String - isotope name 
    + signal::Bool - signal or background
    + activity::Real - activity of the given process in units [Bq]
    + timeMeas::Real - time of measurement in units [s]
    + nTotalSim::Real - total number of originally simulated events from Falaise
    + bins::AbstractRange - binning to be used in format (min:step:max)
    + amount::Real - mass [kg] or volume [l / m³] of the object where the isotope is present: i.e. source foil
    + efficiency::Hist2D - the 2D histogram with efficiencies per bin 
"""
mutable struct Process
    dataVector::Vector{<:Real}
    isotopeName::String
    signal::Bool
    activity::Real
    timeMeas::Real
    nTotalSim::Real
    bins::AbstractRange
    amount::Real
    efficiency::Hist2D
end

"""
    mutable struct Process -> `(dataVector, isotopeNam, signa, activit, timeMea, nTotalSi, bin, efficiency)`

The following fields are defined:
    + dataVector::Vector{<:Real} - vector of initial data points
    + isotopeName::String - isotope name 
    + signal::Bool - signal or background
    + activity::Real - activity of the given process in units [Bq]
    + timeMeas::Real - time of measurement in units [s]
    + nTotalSim::Real - total number of originally simulated events from Falaise
    + bins::AbstractRange - binning to be used in format (min:step:max)
    + amount::Real - mass [kg] or volume [l / m³] of the object where the isotope is present: i.e. source foil
    + efficiency::Hist2D - the 2D histogram with efficiencies per bin 
"""
function Process(
    dataVector::Vector{<:Real},
    isotopeName::String,
    signal::Bool,
    activity::Real,
    timeMeas::Real,
    nTotalSim::Real,
    bins::AbstractRange,
    amount::Real
)
    eff = get_efficiency(dataVector, bins, nTotalSim)
    return Process(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount, eff)
end

function Process(
    dataVector::Vector{<:Real},
    processDict::Dict
)
    @unpack isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount = processDict
    eff = get_efficiency(dataVector, bins, nTotalSim)
    return Process(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount, eff)
end

function get_process(isotopeName, processesVector)
    idx = findfirst( x -> x.isotopeName == isotopeName , processesVector )
    return processesVector[idx]
end

"""
    get_efficiency(dataVector::Vector{<:Real}, bins::AbstractRange, nTotalSim::Real)

    returns a 2D histgoram of type Hist2D with efficiencies per bin.

"""
function get_efficiency(dataVector::Vector{<:Real}, bins::AbstractRange, nTotalSim::Real)
    h2d=get_nPassed(dataVector, bins, wgt=inv(nTotalSim)) # nPassed scaled by total simulated events
    h2d.sumw2 .= sqrt.(h2d.sumw2) # weights scaled 
    return h2d
end

"""
    get_nPassed(process::Process, bins)

    returns a 2D histgoram of type Hist2D with passed events per bin.

"""
function get_nPassed(process::Process, bins::AbstractRange; wgt=1)
    return get_nPassed(process.dataVector, bins, wgt=wgt)
end

function get_nPassed(dataVector::Vector{<:Real}, bins::AbstractRange; wgt=1)
    h1d = Hist1D(dataVector; binedges=bins)
    h2d = Hist2D(; counttype=Float64, binedges=(bins, bins))#prepare empty 2d Histogram

    bc=bincenters(h1d)
    for (i,b) in enumerate(bc)
        cs=cumsum(bincounts(h1d)[i:end]) * wgt
        push!.(h2d, b, bc[i:end], cs)
        # bincounts(h2d)[i, i:end]=cumsum(bincounts(h1d)[i:end])
    end

    h2d
end


"""
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts(processes::Process...)
    h2d = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes
        if (p.signal)
            @warn("get_bkg_counts(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            if (eltype(p.activity) <: Measurement) # check if measurement with uncertainties
                h2d += p.efficiency * p.activity.val * p.amount * p.timeMeas
            else
                h2d += p.efficiency * p.activity * p.amount * p.timeMeas
            end
        end
    end
    return h2d
end

"""
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts_1D(processes::Process...)
    h1d = Hist1D(; binedges=(processes[1].bins))

    for p in processes
        if (p.signal)
            @warn("get_bkg_counts(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            h1d += get_bkg_counts_1D(p)
        end
    end
    return h1d
end

function get_bkg_counts_1D(process::Process)
    fh1 = Hist1D(; binedges=process.bins) 
    push!.(fh1, process.dataVector )

    totalEff = length(process.dataVector) / process.nTotalSim
    fh1 = normalize(fh1, width = false) 
    fh1 *= totalEff * process.amount * process.timeMeas * Measurements.value(process.activity)
    
    return fh1
end

"""
    Returns the number of expected counts of the given process. 
    (n = ε⋅t⋅m⋅A)
"""
function get_sig_counts(processes::Process...)
    h2d = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes
        if (!p.signal)
            @warn("get_sig_counts(): passed isotope $(p.isotopeName) is a background process!!")
        else
            if (eltype(p.activity) <: Measurement) # check if measurement with uncertainties
                h2d += p.efficiency * p.activity.val * p.amount * p.timeMeas
            else
                h2d += p.efficiency * p.activity * p.amount * p.timeMeas
            end
        end
    end
    return h2d
end

"""
    returns epsilon/(S(b)) 

"""
function get_epsilon_to_b(α, processes::Process...; approximate="formula")
    ε = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes # sum efficiencies of the signal processes (does this make sense?)
        if(p.signal)
            ε += p.efficiency
        end
    end

    backgroundCounts = get_bkg_counts(processes...)
    backgroundCounts.bincounts .= get_FC.( backgroundCounts.bincounts , α; approximate=approximate)
    epsToB = ε/(backgroundCounts)

    replace!(epsToB.bincounts, NaN => 0.0)

    return epsToB
end

function get_sToBRatio(processes::Process...)
    backgroundCounts = get_bkg_counts(processes...)
    signalCounts = get_sig_counts(processes...)

    if size(bincounts(backgroundCounts)) != size(bincounts(signalCounts))
        throw(ArgumentError("Input matrices must have the same dimensions"))
    end

    StoB = signalCounts / backgroundCounts
    replace!(StoB.bincounts, NaN => 0.0)

    return StoB
end

function get_tHalf_map(SNparams, α, processes::Process...; approximate="formula")
    ε = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes # sum efficiencies of the signal processes (does this make sense?)
        if(p.signal)
            ε += p.efficiency
        end
    end

    b = get_bkg_counts(processes...)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
     
    b.bincounts .= get_FC.(b.bincounts, α; approximate=approximate)

    tHalf = constantTerm * ε / (b)
    replace!(tHalf.bincounts, NaN => 0.0)

    return tHalf
end

function get_tHalf_map(SNparams, α, processes::Process...; approximate="formula")
    ε = Hist2D(; binedges=(processes[1].bins, processes[1].bins))

    for p in processes # sum efficiencies of the signal processes (does this make sense?)
        if(p.signal)
            ε += p.efficiency
        end
    end

    b = get_bkg_counts(processes...)

    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    constantTerm = log(2) * (Nₐ / W) * (foilMass * a * tYear )
     
    b.bincounts .= get_FC.(b.bincounts, α; approximate=approximate)

    tHalf = constantTerm * ε / (b)
    replace!(tHalf.bincounts, NaN => 0.0)

    return tHalf
end

"""
    get_isotope_details( process::Process )

    returns a tuple of isotope details: (a, m, t, ε, nExpTot, nTotSim)

    Useful for printing and estimating expected counts.
"""
function get_isotope_details( process::Process )
    a = process.activity
    m = process.amount
    t = process.timeMeas
    nTotSim = process.nTotalSim
    ε = length(process.dataVector) / nTotSim

    nExpTot = a*m*t*ε

    return a, m, t, ε, nExpTot, nTotSim 
end

"""
    print_isotope_details( processes::Process... )

    returns a string with basic isotope details.
"""
function print_isotope_details( process::Process )
    (a, m, t, ε, nExpTot, nTotSim) = get_isotope_details(process) 

    a = round(a, sigdigits=3)
    ε = round(ε*100, sigdigits=3)
    m = round(m, sigdigits=3)
    nExpTot = round(nExpTot, sigdigits=3)
    isotope = process.isotopeName
    
    "|$isotope | $nExpTot | $ε | $a | $m |"
end

"""
    get_estimated_bkg_counts(best_ROI, SNparams, processes::Process...)
==========
Outputs the estimated bkg counts for the given ROI. 
This process is obtained from looking up the `bkg_rate` at `best_ROI` and multiplying by mass and time.  
"""
function get_bkg_counts_ROI(minBinCenter, maxBinCenter,processes::Process...)
    b = get_bkg_counts(processes...)

    return lookup(b , minBinCenter, maxBinCenter)
end

function get_bkg_counts_ROI(best_ROI, processes::Process...)
    binStepHalf = step(processes[1].bins) / 2              # get the binning step and divide by half
    minBinCenter = best_ROI[:minBinEdge] + binStepHalf       # get the center of the minimal bin in ROI
    maxBinCenter = best_ROI[:maxBinEdge] - binStepHalf       # get the center of the maximal bin in ROI

    return get_bkg_counts_ROI(minBinCenter, maxBinCenter, processes...)
end

function FHist.lookup(process::Process, best_ROI::Dict) 
    binStepHalf = step(process.bins) / 2              # get the binning step and divide by half
    minBinCenter = best_ROI[:minBinEdge] + binStepHalf       # get the center of the minimal bin in ROI
    maxBinCenter = best_ROI[:maxBinEdge] - binStepHalf       # get the center of the maximal bin in ROI
    
    return lookup(process.efficiency, minBinCenter, maxBinCenter)
end

FHist.lookup(process::Process, x::Real, y::Real) = lookup(process.efficiency, x, y)
DrWatson.default_allowed(::Process) = (Real, String, Bool, AbstractRange)
DrWatson.allaccess(::Process) = (:isotopeName, :signal, :bins, :activity, :nTotalSim, :amount)



