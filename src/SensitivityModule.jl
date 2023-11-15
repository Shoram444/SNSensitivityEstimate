module SensitivityModule

using StatsPlots, DataFramesMeta, LaTeXStrings, FHist, StatsBase, UnROOT, DrWatson


include("Misc.jl")
export 
    generate_raw_plots,
    fill_from_root_file,
    get_max_bin

include("Process.jl")
export 
    Process,
    get_efficiency,
    get_nPassed,
    get_bkg_rate,           # obtain the rate of bkg processes per ROI: b_rate = sum_i A_i * eff_i  (where i is the individual bkg process)
    get_sig_rate,           # obtain the rate of sig processes per ROI: s_rate = sum_i A_i * eff_i  (where i is the individual sig process)
    get_sToBRatio           # sig_rate / bkg_rate



end # module
