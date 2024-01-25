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
            if( eltype(p.activity) <: Measurement ) # check if measurement with uncertainties
                h2d += p.efficiency * p.activity.val * p.amount * p.timeMeas
            else
                h2d += p.efficiency * p.activity * p.amount * p.timeMeas
            end
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
            if( eltype(p.activity) <: Measurement ) # check if measurement with uncertainties
                h2d += p.efficiency * p.activity.val * p.amount * p.timeMeas
            else
                h2d += p.efficiency * p.activity * p.amount * p.timeMeas
            end
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

"""
    get_estimated_bkg_counts(best_ROI, SNparams, processes::Process...)
==========
Outputs the estimated bkg counts for the given ROI. 
This process is obtained from looking up the `bkg_rate` at `best_ROI` and multiplying by mass and time.  
"""
function get_estimated_bkg_counts(minBinCenter, maxBinCenter, SNparams, processes::Process...)
    for p in processes
        if(p.signal)
            error("Isotope $(p.isotopeName) is a signal process! Please provide only background processes!")
        end
    end
    h2d = Hist2D(Float64, bins=(processes[1].bins, processes[1].bins))
    bkg_cts_mat = zeros(size(bincounts(h2d)))   

    for p in processes
        if (p.signal)
            @warn("get_bkg_rate(): passed isotope $(p.isotopeName) is a signal process!!")
        else
            bkg_cts_mat += bincounts(p.efficiency) * p.activity * p.amount * p.timeMeas
        end
    end

    return lookup(bkg_cts_mat, minBinCenter, maxBinCenter, processes[1].bins)
end

function get_estimated_bkg_counts(best_ROI, SNparams, processes::Process...)
    binStepHalf = step( processes[1].bins ) / 2              # get the binning step and divide by half
    minBinCenter = best_ROI[:minBinEdge] + binStepHalf       # get the center of the minimal bin in ROI
    maxBinCenter = best_ROI[:maxBinEdge] - binStepHalf       # get the center of the maximal bin in ROI

    return get_estimated_bkg_counts(minBinCenter, maxBinCenter, SNparams, processes...)
end


DrWatson.default_allowed(::Process) = (Real, String, Bool, AbstractRange)
DrWatson.allaccess(::Process) = (:isotopeName, :signal, :bins, :activity, :nTotalSim, :amount)

function FHist.lookup(weights::Matrix, x::Real, y::Real, binning)
    h2d = Hist2D(Float64, bins=(binning, binning))
    h2d.hist.weights = [bc.val for bc in weights]

    rx, ry = binedges(h2d)
    !(first(rx) <= x <= last(rx)) && return missing
    !(first(ry) <= y <= last(ry)) && return missing

    return weights[FHist._edge_binindex(rx, x), FHist._edge_binindex(ry, y)]
end

