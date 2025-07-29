using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, FHist, ColorSchemes


include(scriptsdir("ND/results/best_rois.jl"))
# include(srcdir("params/Params.jl"))

analysisDict    = Dict(:signal => "bb0nu_foil_bulk", :bining => (300, 3500), :bin_width => 100, :mode => "sumE")

roi             = bb0nu_roi
roi[:sumE]      = (2700, 3100)


Bin_low, Bin_high, bin_width = analysisDict[:bining][1], analysisDict[:bining][2], analysisDict[:bin_width]
backgrounds     = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_wire_surface", "Tl208_foil_bulk", "gamma_experimental_surface"]

all_processes   = load_data_processes(datadir("mva/fal5_TKrec_J40"), analysisDict[:mode], fwhm=0.0, roi = roi)
signal          = get_process(analysisDict[:signal], all_processes) |> first
background      = [get_process(b, all_processes) |> first for b in backgrounds]


set_nTotalSim!(signal, 0.1e8)
set_nTotalSim!(background[1], 0.1e8)
set_activity!(background[1], 0.0019235366786346892)
set_activity!(background[end], 0.185)


for b in background
    set_bins!(b, Bin_low:bin_width:Bin_high+bin_width)
end
set_bins!(signal, Bin_low:bin_width:Bin_high+bin_width)

set_signal!(background[1], false)

bkgs = [get_bkg_counts_1D(b) for b in background]
bkg_sums = round.(integral.(bkgs), digits = 3)
activities = [b.activity for b in background]
bkg_errs = [round(sum(sqrt.(sumw2(b))), digits=3) for b in get_bkg_counts_1D.(background)]

using PrettyTables, DataFramesMeta, Measurements

pretty_table(
    DataFrame(
        process = [L"$2\nu\beta\beta$", L"$^{214}$Bi",L"Radon $$", L"$^{208}$Tl", L"external $$", "total"],
        exp_bkg = vcat(bkg_sums .± bkg_errs, [sum(bkg_sums .± bkg_errs)]),
        activities = vcat(activities, [0.0]),
    ),
    header = ["source", "counts in ROI", "activitiies"],
    backend = Val(:markdown)
)
