"""
    halfLife_to_activity( NA::Real, W::Real, Thalf::Real ) -> returns activity in [Bq/kg]
"""
function halfLife_to_activity(NA::Real, W::Real, Thalf::Real)
    return (log(2) * NA) / (W * Thalf)
end

"""
    get_tHalf(SNparams, efficiency, bkgCounts, α; approximate = :true) -> returns the sensitivity in yr. 
"""
function get_tHalf(SNparams, efficiency, bkgCounts, α; approximate=:true)
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    if approximate
        b = α√bkgCounts
    end

    tHalf = log(2) * (Nₐ / W) * efficiency * (foilMass * a * tYear / b)
end


"""
    get_isotope_details( activityParams, SNParams, simParams, isotope, dfData )

    returns a tuple of isotope details: (a, m, t, ε, nExpTot, nTotSim)

    Useful for printing and estimating expected counts.
"""
function get_isotope_details( activityParams, SNParams, simParams, isotope, dfData )
    a = activityParams[isotope]
    m = if( occursin("foil_bulk", string(isotope)))
        SNParams["foilMass"]
    elseif( occursin("surface", string(isotope)))
        SNParams["gasVolume"]
    elseif( occursin("PMT_bulk", string(isotope)))
        SNParams["PMTGlassMass"]
    elseif( occursin("wire_bulk", string(isotope)))
        SNParams["wireBulkMass"]
    end
    t = SNParams["t"]
    nTotSim = simParams[isotope]
    ε = nrow(dfData) / nTotSim

    nExpTot = a*m*t*ε

    return a, m, t, ε, nExpTot, nTotSim 
end

"""
    print_isotope_details( activityParams, SNParams, simParams, isotope, dfData )

    returns a string with basic isotope details.
"""
function print_isotope_details( activityParams, SNParams, simParams, isotope, dfData )
    (a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( activityParams, SNParams, simParams, isotope, dfData ) 

    a = round(a, sigdigits=3)
    ε = round(ε*100, sigdigits=3)
    m = round(m, sigdigits=3)
    nExpTot = round(nExpTot, sigdigits=3)
    
    "|$isotope | $nExpTot | $ε | $a | $m |"
end

"""
    estimated_counts_hist1D(dataDf, activity, amount, time, binning, nTotalSim)

    returns `Hist1D` histogram with estimated counts of the given process. 
"""
function estimated_counts_hist1D(dataDf, activity, amount, time, binning, nTotalSim)
    totEff = nrow(dataDf) / nTotalSim
    h1d = 
        Hist1D(
            dataDf.reconstructedEnergy1 .+ dataDf.reconstructedEnergy2 ,
            binning
        )
    
    h1d = normalize(h1d; width =:true) # normalize to area of 1
    
    if(typeof(activity)<:Measurement) # just a precaution, sometime I use activities with uncertainty
        activity = activity.val
    end
    
    h1d.hist.weights = h1d.hist.weights .*(
        activity * time * amount * totEff    
    ) # multiply each bin by activity * time * amount * totalEfficiency

    return h1d
end