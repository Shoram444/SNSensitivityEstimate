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
    
    savefig(hPhi, joinpath(plotsdir(), "$(isotope)_raw_angdist.png") )
    savefig(hEne, joinpath(plotsdir(), "$(isotope)_raw_enedist.png") )

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