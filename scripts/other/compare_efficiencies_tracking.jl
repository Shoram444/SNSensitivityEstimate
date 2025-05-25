using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, Revise, FHist, PrettyTables, DataFramesMeta, LinearAlgebra


# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(srcdir("params/Params.jl"))

# Dictionary with the analysis parameters. (Nice for when you want to save stuff and make sure you know what params you used in this analysis.)
analysisDict = Dict(
    :Bfield => "Boff", # magnetic field on/off
    :Eres => "8perc", # FWHM of the OMs (sorry for the naming...)
    :mode => "singleE", 
    :trackAlgo => "TKrec",
    :signal => "bb0nu_foil_bulk", # signal process
    :neutron_config => "no_neutron"
)

# files_directory = "fal5_$(analysisDict[:Eres])_$(analysisDict[:Bfield])_$(analysisDict[:trackAlgo])_twoDistinct_edep_bcu"
# files_directory = "fal5_8perc_Boff_CAT_evis_bcu_J39"
files_directory = "fal5_8perc_Boff_TKrec_evis_bcu_J38"



# Load all the processes in the directory. Function `load_processes` takes two arguments:
#  1. dir::String -> the name of the directory where the root files are stored
#  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
data_processes = load_data_processes(
    files_directory, 
    analysisDict[:mode],
    fwhm = 0.0
)

signal = get_process("bb0nuM2_foil_bulk", data_processes) |> first
# signal = get_process("RH037_foil_bulk", data_processes) |> first

# declare background processes
background = [
    get_process("bb_foil_bulk", data_processes) |> first,
    get_process("Bi214_foil_bulk", data_processes) |> first,
    get_process("Bi214_wire_surface", data_processes) |> first,
    get_process("Tl208_foil_bulk", data_processes) |> first,
    get_process("K40_foil_bulk", data_processes) |> first,
    get_process("Pa234m_foil_bulk", data_processes) |> first,
]


set_signal!(background[1], false)

set_nTotalSim!( signal, 2*1e8 )
set_nTotalSim!( background[1], 2*1e8 )
set_nTotalSim!( background[2], 2*1e8 )
set_nTotalSim!( background[3], 2*1e8 )
set_nTotalSim!( background[4], 2*1e8 )
set_nTotalSim!( background[5], 2*1e8 )
set_nTotalSim!( background[6], 2*1e8 )

println("Processes initialized.")

α = 1.64485362695147

t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
@show best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
@show ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α; approximate = "table"), sigdigits=3)

bkg_hists = get_bkg_counts_1D.(background)
sig_hist = get_bkg_counts_1D(signal)



# labels = [
#     "0nubb",
#     "2nubb",
#     "Bi214",
#     "Radon",
#     "Tl208",
#     "K40",
#     "Pa234m"
# ]

# effs_ROI = [
#     lookup(signal, best_t12ESum),
#     lookup(background[1], best_t12ESum),
#     lookup(background[2], best_t12ESum),
#     lookup(background[3], best_t12ESum),
#     lookup(background[4], best_t12ESum),
#     lookup(background[5], best_t12ESum),
#     lookup(background[6], best_t12ESum)
# ] .* 100

# full_ROI = Dict(
#     :minBinEdge => 0.0,
#     :maxBinEdge => 3500.0,
# )

# effs_full = [
#     lookup(signal, full_ROI),
#     lookup(background[1], full_ROI),
#     lookup(background[2], full_ROI),
#     lookup(background[3], full_ROI),
#     lookup(background[4], full_ROI),
#     lookup(background[5], full_ROI),
#     lookup(background[6], full_ROI)
# ] .* 100

# df = DataFrame(
#     process = labels,
#     efficiency_ROI = effs_ROI,
#     efficiency_full = effs_full,
# )

# pt = pretty_table(
#     df,
#     header = (["Process", "Efficiency (ROI)", "Efficiency (Full)"], [" ", "%", "%"]),
#     backend = Val(:markdown),
#     formatters = ft_printf("%.3E"),
# )

# sn = savename("eff", analysisDict, "md")

# open(scriptsdir("other","efficiencies",sn), "w") do io
#     pretty_table(
#         io, 
#         df,
#         header = (["Process", "Efficiency (ROI)", "Efficiency (Full)"], [" ", "%", "%"]),
#         backend = Val(:markdown),
#         formatters = ft_printf("%.3E"),
#     )
# end