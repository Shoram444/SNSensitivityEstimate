using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta


# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
include(scriptsdir("Params.jl"))

# Dictionary with the analysis parameters. (Nice for when you want to save stuff and make sure you know what params you used in this analysis.)
analysisDict = Dict(
    :Bfield => "Boff", # magnetic field on/off
    :Eres => "8perc", # FWHM of the OMs (sorry for the naming...)
    :mode => "sumE", 
    :trackAlgo => "TIT",
    :signal => "bb0nu",
)


# Load all the processes in the directory. Function `load_processes` takes two arguments:
#  1. dir::String -> the name of the directory where the root files are stored
#  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
all_processes = load_processes(
    "fal5_$(analysisDict[:Eres])_$(analysisDict[:Bfield])_$(analysisDict[:trackAlgo])_twoDistinct", 
    analysisDict[:mode]
)

# declare which process is signal
signal = get_process("bb0nu_foil_bulk", all_processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 1e8 )
set_nTotalSim!( background[2], 1e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 1e8 )

println("Processes initialized.")

# For now I only implemented 90% CL estimation, maybe someday I'll get to more...
α = 1.64485362695147

# The base structure of the `SensitivityModule` is the `Process`. For details, take a look into `src/Process.jl`, or hit `? Process` in REPL.
# The way this works is that by loading a process, we create an efficiency map for each isotope (process). For description of what these maps represent
# see: #docdb-5833. Stop at page 8, because there was a mistake in the presentation. 
# (The mistake: it's not enough to calculate efficiency ratios, you need to calculate sensitivity itself...)
# - The efficiency maps for all background processes in the set are combined into a full `expected background counts map`.
# - Together with the signal efficiency map, the function `get_tHalf_map` calculates the sensitivity for each ROI in the map.
# 
# Having the map (which is essentially a 2D histogram), we pick the highest value (and the bin edges) with `get_max_bin`.
# That's basically it, if you're only interested in sensitvity.
t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
best_t12ESum = get_max_bin(t12MapESum)

# If you want additional info, like background counts in the ROI, use: `get_bkg_counts_ROI`.
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)

# To get the signal efficiency at the ROI, use: `lookup(signal, best_t12ESum)`
effbb = lookup(signal, best_t12ESum)


# Plotting:

# 2D tHalf map:
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
      b  = $(round(expBkgESum, sigdigits = 3)) 
      T12 ≥  $(round(ThalfbbESum, sigdigits=2)) yr 
      ε = $(round(effbb, sigdigits = 3)*100)%"

let 
    f = Figure(size=(600, 400))
    a = Axis(f[1,1], xlabel = "min_ROI [keV]", ylabel = "max_ROI [keV]")
    p = plot!(a, t12MapESum)
    text!(a, 2000, 500, text=lbl)
    Colorbar(f[1,2], p, label="sensitivity [yr]", scale=log10)
    saveName = savename("tHalf_map", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f
end


# Sensitivity as a function of detector life-time:
t = range(0, 6, 60)

t12(t, e, b) = get_tHalf(
    SNparams["W"],
    SNparams["foilMass"],
    SNparams["Nₐ"],
    t,
    SNparams["a"],
    e,
    b / SNparams["tYear"] * t,
    α;
    approximate="formula"
)

let 
    f = Figure(size=(600, 400))
    a = Axis(f[1,1], xlabel = "detector life-time [yr]", ylabel = "sensitivity [yr]", limits= (0,6, nothing, nothing))
    p = lines!(a, t, t12.(t, effbb,expBkgESum), label = "FC: 8% resolution; simulation")
    axislegend(a, position = :lt)
    saveName = savename("sensitivity_in_time", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f
end


# Total background model:
bkg_hists = get_bkg_counts_1D.(background)

with_theme(theme_latexfonts()) do
    f = Figure()
    ax = Axis(f[1,1], xlabel = analysisDict[:mode], ylabel = "counts", yscale = log10, limits = (nothing, nothing, 1e-5, 1e6))
    
    colors = colorschemes[:tableau_colorblind]
    labels = [b.isotopeName for b in background]
    st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
	for i=2:length(bkg_hists)
		hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
	end
    
    ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
    axislegend(ax, position = :lt, orientation = :horizontal, nbanks=2)
    saveName = savename("background_model", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f
end

# Background by source:
ROI_a, ROI_b = best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]

bkgs = [sum(bincounts(restrict(b, ROI_a, ROI_b)))  for b in get_bkg_counts_1D.(background)]

pretty_table(
    DataFrame(
        process = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_radon", "Tl208_foil_bulk", "total"],
        counts = vcat(bkgs, sum(bkgs)),
        activity_used = vcat([background[i].activity for i in 1:4], "--"),
    ),
    backend = Val(:markdown),
)
