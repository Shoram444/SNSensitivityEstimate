function generate_raw_plots(inDf::DataFrame, isotope; kwargs...)
    phi, e1, e2 = inDf[!, :phi], inDf[!, :reconstructedEnergy1], inDf[!, :reconstructedEnergy2]

    plot_args = Dict(
        :left_margin => 12Plots.mm,
        :right_margin => 12Plots.mm,
        :top_margin => 12Plots.mm,
        :bottom_margin => 12Plots.mm,
        :size => (1200, 900),
        :dpi => 200,
        :titlefontsize => 24,
        :guidefontsize => 22,
        :tickfontsize => 20,
        :legendfontsize => 24,
        :formatter => :auto,
        :legend => :best,
        :lw => 4,
        :thickness_scaling => 1.1
    )

    hPhi = stephist(
        phi;
        nbins=0:5:180,
        xlabel="escape angle "*L"[\degree]",
        ylabel="counts/" * L"5\degree",
        title= "angular distribution of $isotope",
        label = "",
        plot_args...
    )

    hEne = stephist(
        e1;
        nbins=0:100:3500,
        xlabel="single-electron energy [keV]",
        ylabel="counts/" * "100 keV",
        title= "single electron energy distribution of $isotope",
        label = "energy 1",
        plot_args...
        )

    stephist!(
        e2;
        nbins=0:100:3500,
        label = "energy 2",
        plot_args...
    )

    stephist!(
        vcat(e1,e2);
        nbins=0:100:3500,
        label = "both electrons",
        plot_args...
    )
    
    safesave(plotsdir( "Raw", "$(isotope)_raw_angdist.png"), hPhi )
    safesave(plotsdir("Raw", "$(isotope)_raw_enedist.png"), hEne)

    plot(hPhi, hEne, size = (1800, 600))

end

function fill_from_root_file(inFile::ROOTFile, treeName::String, fieldNames)
    df = LazyTree(inFile, treeName, [fieldNames...]) |> DataFrame

    for k in fieldNames
        filter!(string(k) => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df)
    end
    return df
end

FHist.binedges(h::Histogram) = h.edges[1]


"""
    get_max_bin(h2d::Hist2D)

Returns a Dict with the following keys:
    + :minBinEdge => minBinCenterEdge - binStepHalf ->  defines the lower edge of ROI
    + :maxBinEdge => maxBinCenterEdge + binStepHalf ->  defines the upper edge of ROI
    + :maxBinCount => maxBinCount ->  provides the maximum bin counts of the 2D histogram
"""
function get_max_bin(h2d::Hist2D)
    (step(binedges(h2d)[1]) != step(binedges(h2d)[1])) && error("bins must be the same!")
    binStepHalf = step(binedges(h2d)[1])/2
    BinID = argmax(bincounts(h2d))
    minBinEdge, maxBinEdge = bincenters(h2d)[1][BinID[1]], bincenters(h2d)[1][BinID[2]]
    maxBinCount = lookup(h2d, minBinEdge, maxBinEdge)

    return Dict(
        :minBinEdge => minBinEdge - binStepHalf,
        :maxBinEdge => maxBinEdge + binStepHalf,
        :maxBinCount => maxBinCount
    )
end


"""
    halfLife_to_activity( NA::Real, W::Real, Thalf::Real ) -> returns activity in [Bq/kg]
"""
function halfLife_to_activity( NA::Real, W::Real, Thalf::Real )
    return (log(2) * NA) / (W * Thalf) 
end

