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
    DensityInterface,
    Combinatorics

# Define main module abstract type


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

include("AbstractProcess.jl")
export  
    AbstractProcess,
    get_process,
    set_signal!,
    get_bkg_counts,
    get_bkg_counts_1D,
    get_sig_counts,
    get_epsilon_to_b,
    get_sToBRatio,
    get_tHalf_map,
    get_bkg_counts_ROI

include("HistProcess.jl")
export 
    HistProcess,
    get_counts,
    get_nPassed,
    get_bkg_counts,
    get_bkg_counts_1D,
    get_isotope_details,
    print_isotope_details,
    get_process

include("DataProcess.jl")
export 
    DataProcess,
    set_activity!,
    set_timeMeas!,
    set_nTotalSim!,
    set_bins!,
    set_amount!,
    get_efficiency,
    get_nPassed,
    get_bkg_counts_1D,
    get_bkg_counts
    
include("DataProcess3D.jl")
export
    SensitivityEstimate3D,
    ROIEfficiency3D,
    DataProcess3D,
    get_effciency3D,
    get_bkg_counts,
    get_sensitivity3D,
    get_s_to_b

include("NDimDataProcess.jl")
    export
        SensitivityEstimateND,
        ROIEfficiencyND,
        DataProcessND,
        get_roi_effciencyND,
        get_roi_bkg_counts,
        get_s_to_b,
        get_sensitivityND

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
    load_data_processes,
    load_hist_processes,
    load_3D_processes#,
    # load_ndim_processes

include("Models.jl")
export 
    turing_normal_exponential_model,
    BAT_normal_exponential_model 
    
end # module
