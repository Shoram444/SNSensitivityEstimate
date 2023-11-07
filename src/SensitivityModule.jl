module SensitivityModule

using StatsPlots, DataFramesMeta, LaTeXStrings, FHist, StatsBase, UnROOT, DrWatson


include("Misc.jl")
export 
    generate_raw_plots,
    fill_from_root_file

include("Process.jl")
export 
    Process,
    plot_efficiency,
    get_efficiency,
    get_nPassed,
    plot_nPassed,
    lookup_efficiency,
    get_bkg_counts,
    get_sig_counts,
    get_sToBRatio,
    plot_sToBRatio



end # module
