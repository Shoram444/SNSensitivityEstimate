using CSV, DataFrames

"""
    halfLife_to_activity( NA::Real, W::Real, Thalf::Real ) -> returns activity in [Bq/kg]
"""
function halfLife_to_activity(NA::Real, W::Real, Thalf::Real)
    return (log(2) * NA) / (W * Thalf)
end

"""
    get_tHalf(SNparams, efficiency, bkgCounts, α; approximate = :true) -> returns the sensitivity in yr. 
"""
function get_tHalf(SNparams, efficiency, bkgCounts, α; approximate="formula")
    @unpack W, foilMass, Nₐ, tYear, a = SNparams
    b = get_FC(bkgCounts, α; approximate=approximate)

    tHalf = log(2) * (Nₐ / W) * (foilMass * a * tYear) * efficiency / b
end

function get_tHalf(W, foilMass, Nₐ, tYear, a , efficiency, bkgCounts, α; approximate="formula")
    b = get_FC(bkgCounts, α; approximate=approximate)
    tHalf = log(2) * (Nₐ / W) * (foilMass * a * tYear) * efficiency / b
end

function get_FC(b, α; approximate="table")
    b = Measurements.value(b) # need to figure out how to work this with BBO
    if( approximate == "formula")
        if (b < 30 && isapprox(1.64, α, atol=0.1))
            b = 2.44+0.8467*b-0.08849*b^2+0.00625*b^3-0.0002124*b^4+0.000002712*b^5
        else
            b = α * sqrt(b)
        end
    elseif( approximate == "table")
        if (b < 329 && isapprox(1.64, α, atol=0.1))
            tbl = CSV.read(srcdir("MPFC_table90.csv"), DataFrame)
            b = tbl[findfirst(tbl[:,1].== round(b, digits =2)), 2]
        else 
            b = α * sqrt(b)
        end
    else
        b = α * sqrt(b)
    end
    return b
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
    elseif( occursin("PMT_glass_bulk", string(isotope)))
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
            dataDf.reconstructedEnergy1 .+ dataDf.reconstructedEnergy2;
            binedges = binning
        )
    
    h1d = normalize(h1d; width =:false) # normalize to area of 1
    
    if(typeof(activity)<:Measurement) # just a precaution, sometime I use activities with uncertainty
        activity = activity.val
    end
    
    h1d.bincounts .= h1d.bincounts .*(
        activity * time * amount * totEff    
    ) # multiply each bin by activity * time * amount * totalEfficiency

    return h1d
end