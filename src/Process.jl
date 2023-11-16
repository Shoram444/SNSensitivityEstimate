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
    + efficiency::Hist2D - the 2D histogram with efficiencies per bin 
"""
function Process(
    dataVector::Vector{<:Real}, 
    isotopeName::String, 
    signal::Bool, 
    activity::Real, 
    timeMeas::Real, 
    nTotalSim::Real, 
    bins::AbstractRange
) where {T<:Real}
    eff = get_efficiency(dataVector, bins, nTotalSim)
    return Process(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins, eff)
end

"""
    get_efficiency(dataVector::Vector{<:Real}, bins::AbstractRange, nTotalSim::Real)

    returns a 2D histgoram of type Hist2D with efficiencies per bin.

"""
function get_efficiency(dataVector::Vector{<:Real}, bins::AbstractRange, nTotalSim::Real)
    n = nTotalSim

    h2d = get_nPassed(dataVector, bins) * (1 / n) # Matrix - scalar multiplication
end

"""
    get_nPassed(process::Process, bins)

    returns a 2D histgoram of type Hist2D with passed events per bin.

"""
function get_nPassed(process::Process, bins::AbstractRange)
    return get_nPassed(process.dataVector, bins)
end

function get_nPassed(dataVector::Vector{<:Real}, bins::AbstractRange)
    h1d = Hist1D(dataVector, bins)
    h2d = Hist2D(Float64, bins=(bins, bins))#prepare empty 2d Histogram

    binCenters = collect(bincenters(h1d))
    binStep = step(bins)
    for (i, xVal) in enumerate(binCenters)
        # for (j, yVal) in enumerate(binCenters[i+1:end])
        for (j, yVal) in enumerate(binCenters[i:end])
            # j += i
            integral = sum(lookup.(h1d, xVal:binStep:yVal))
            push!(h2d, xVal, yVal, integral)
        end
    end

    h2d
end


FHist.lookup(process::Process, x::Real, y::Real) = lookup(process.efficiency, x, y)

function get_bkg_rate(processes::Process...)
    h2d = Hist2D(Float64, bins=(processes[1].bins, processes[1].bins))

    for p in processes
        if (p.signal)
            @warn("get_bkg_rate(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            h2d += p.efficiency * p.activity
        end
    end
    return h2d
end

function get_sig_rate(processes::Process...)
    h2d = Hist2D(Float64, bins=(processes[1].bins, processes[1].bins))

    for p in processes
        if (!p.signal)
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

DrWatson.default_allowed(::Process) = (Real, String, Bool, AbstractRange)
DrWatson.allaccess(::Process) = (:isotopeName, :signal, :bins)
