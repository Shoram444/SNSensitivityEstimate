module SensitivityModule

using StatsPlots, DataFramesMeta, LaTeXStrings, FHist, StatsBase, UnROOT, DrWatson, MakieCore


include("Misc.jl")
export 
    generate_raw_plots,
    fill_from_root_file,
    get_max_bin,
    halfLife_to_activity,
    get_tHalf,
    annotatewithbox!,
    get_2D_vertex_separation,
    get_1D_vertex_separation,
    add_vertex_2D_separation_column!,
    add_vertex_dy_separation_column!,
    add_vertex_dz_separation_column!

include("Process.jl")
export 
    Process,
    get_efficiency,
    get_nPassed,
    get_bkg_rate,               # obtain the rate of bkg processes per ROI: b_rate = sum_i A_i * eff_i  (where i is the individual bkg process)
    get_sig_rate,               # obtain the rate of sig processes per ROI: s_rate = sum_i A_i * eff_i  (where i is the individual sig process)
    get_sToBRatio,              # sig_rate / bkg_rate
    get_estimated_bkg_counts    # obtain the estimated bkg numbers (scaled to mass and Tmeas)

include("PlotsRecipes.jl")

include("Makie_themes.jl")
export 
    dao 



end # module
