

"""
    struct DataProcess -> `(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins, efficiency)`

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

mutable struct DataProcess <: AbstractProcess
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
    mutable struct DataProcess -> `(dataVector, isotopeNam, signal, activity, timeMeas, nTotalSim, bins, efficiency)`

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
function DataProcess(
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
    return DataProcess(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount, eff)
end

function DataProcess(
    dataVector::Vector{<:Real},
    processDict::Dict
)
    @unpack isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount = processDict
    eff = get_efficiency(dataVector, bins, nTotalSim)
    return DataProcess(dataVector, isotopeName, signal, activity, timeMeas, nTotalSim, bins, amount, eff)
end



############################################################
########### setters!
############################################################


function set_activity!(process::DataProcess, activity::Real)
    process.activity = activity
    return process
end

function set_timeMeas!(process::DataProcess, timeMeas::Real)
    process.timeMeas = timeMeas
    return process
end

function set_nTotalSim!(process::DataProcess, nTotalSim::Real)
    process.nTotalSim = nTotalSim
    process.efficiency = get_efficiency(process.dataVector, process.bins, nTotalSim)
    return process
end

function set_bins!(process::DataProcess, bins)
    process.bins = bins
    process.efficiency = get_efficiency(process.dataVector, bins, process.nTotalSim::Real)
    return process
end


function set_amount!(process::DataProcess, amount::Real)
    process.amount = amount
    return process
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
!!
    get_nPassed(process::DataProcess, bins)

    returns a 2D histgoram of type Hist2D with passed events per bin.

"""
function get_nPassed(process::DataProcess, bins::AbstractRange; wgt=1)
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
!! - here it's free for HistProcess
    Returns the number of expected counts of the given process. 
    ( n = (ε⋅t⋅m⋅A) ); 
"""
function get_bkg_counts_1D(process::DataProcess)
    fh1 = Hist1D(; binedges=process.bins) 
    push!.(fh1, process.dataVector )

    totalEff = length(process.dataVector) / process.nTotalSim
    fh1 = normalize(fh1, width = false) 
    fh1 *= totalEff * process.amount * process.timeMeas * Measurements.value(process.activity)
    
    return fh1
end


function get_bkg_counts(process::DataProcess)
    h2d = Hist2D(; binedges=(process.bins, process.bins))

    if (eltype(process.activity) <: Measurement) # check if measurement with uncertainties
        merge!(h2d, process.efficiency * process.activity.val * process.amount * process.timeMeas)
    else
        merge!(h2d, process.efficiency * process.activity * process.amount * process.timeMeas)
    end
    return h2d
end


function get_sig_counts(process::DataProcess)
    h2d = Hist2D(; binedges=(process.bins, process.bins))
    if (eltype(process.activity) <: Measurement) # check if measurement with uncertainties
        merge!(h2d, process.efficiency * process.activity.val * process.amount * process.timeMeas)
    else
        merge!(h2d, process.efficiency * process.activity * process.amount * process.timeMeas)
    end
    return h2d
end




