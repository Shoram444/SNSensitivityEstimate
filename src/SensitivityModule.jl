module SensitivityModule

using 
    StatsPlots, 
    DataFramesMeta, 
    LaTeXStrings, 
    FHist, 
    StatsBase, 
    UnROOT, 
    DrWatson, 
    MakieCore, 
    Measurements,
    Turing,
    CairoMakie,
    BAT,
    DensityInterface


include("Misc.jl")
export 
    generate_raw_plots,
    fill_from_root_file,
    get_2D_vertex_separation,
    get_1D_vertex_separation,
    add_vertex_2D_separation_column!,
    add_vertex_dy_separation_column!,
    add_vertex_dz_separation_column!,
    generate_pseudo_data,
    get_sigma_keV,
    get_sigma_MeV

include("ExpectedCounts.jl")
export 
    halfLife_to_activity,
    get_tHalf,
    get_isotope_details,
    print_isotope_details,
    estimated_counts_hist1D,
    get_FC

include("HistFuncs.jl")
export 
    get_max_bin

include("Process.jl")
export 
    Process,
    get_process,                # returns the Process by isotopeName from vector of processes
    get_efficiency,
    set_signal!,                # experimental feature!
    set_activity!,
    set_timeMeas!,
    set_nTotalSim!,
    set_bins!,
    set_amount!,
    get_nPassed,
    get_bkg_counts,               # obtain the expected number of bkg processes per ROI: b_rate = α * sqrt( sum_i A_i * eff_i * t * amount) (where i is the individual bkg process)
    get_sig_counts,               # obtain the rate of sig processes per ROI: s_rate = sum_i A_i * eff_i  (where i is the individual sig process)
    get_sToBRatio,              # sig / bkg
    get_epsilon_to_b,           # epsilon / bkg
    get_tHalf_map,              # returns the Hist2D object with T1/2 per ROI
    get_bkg_counts_ROI,          # obtain the estimated bkg numbers in given ROI (scaled to mass and Tmeas)
    get_bkg_counts_1D,             # returns Hist1D of expected backgrounds
    get_isotope_details,
    print_isotope_details
    
include("PlotsRecipes.jl")

include("Makie.jl")
export 
    dao 

include("BayessUtils.jl")
export 
    get_mean_conf,
    get_median_conf   

include("LoadData.jl")
export 
    load_files,
    load_processes

include("Models.jl")
export 
    exponential_normal_model,           
    exponential_uniform_normal_model,
    mixture_model,
    exponential_normal_mixture,
    BAT_normal_exponential_model 
    
end # module
