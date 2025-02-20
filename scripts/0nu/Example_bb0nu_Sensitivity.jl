using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta
using Revise
Revise.track(SensitivityModule)

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
    :signal => "bb0nuM1",
    :neutron_config => "full_shielding"
)

files_directory = "fal5_$(analysisDict[:Eres])_$(analysisDict[:Bfield])_$(analysisDict[:trackAlgo])_twoDistinct_edepbcu"

# Load all the processes in the directory. Function `load_processes` takes two arguments:
#  1. dir::String -> the name of the directory where the root files are stored
#  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
data_processes = load_data_processes(
    files_directory, 
    analysisDict[:mode]
)

hist_processes = load_hist_processes(
    files_directory,  
    analysisDict[:mode]
)


# declare which process is signal
# signal = get_process("bb0nu_foil_bulk", data_processes)
# signal = get_process("bb0nuM1_foil_bulk", data_processes)
signal = get_process("bb0nuM2_foil_bulk", data_processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", data_processes),
    get_process("Bi214_foil_bulk", data_processes),
    get_process("Bi214_wire_surface", data_processes),
    get_process("Tl208_foil_bulk", data_processes),
    get_process("K40_foil_bulk", data_processes),
    get_process("Pa234m_foil_bulk", data_processes),
    get_process("neutron_external", hist_processes, "full shielding no floor flux")
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
# set_nTotalSim!( signal, 0.98e8 )
set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 0.99e8 )
set_nTotalSim!( background[2], 0.96e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 0.76e8 )
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )

println("Processes initialized.")

# For now I only implemented 90% CL estimation, maybe someday I'll get to more...
α = 1.64485362695147

t = range(0, 5, 100)


# The base structure of the `SensitivityModule` is the `Process`. For details, take a look into `src/Process.jl`, or hit `? Process` in REPL.
# The way this works is that by loading a process, we create an efficiency map for each isotope (process). For description of what these maps represent
# see: #docdb-5833. Stop at page 8, because there was a mistake in the presentation. 
# (The mistake: it's not enough to calculate efficiency ratios, you need to calculate sensitivity itself...)
# - The efficiency maps for all background processes in the set are combined into a full `expected background counts map`.
# - Together with the signal efficiency map, the function `get_tHalf_map` calculates the sensitivity for each ROI in the map.
# 
# Having the map (which is essentially a 2D histogram), we pick the highest value (and the bin edges) with `get_max_bin`.
# That's basically it, if you're only interested in sensitvity.
t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
best_t12ESum = get_max_bin(t12MapESum)

# If you want additional info, like background counts in the ROI, use: `get_bkg_counts_ROI`.
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)

# To get the signal efficiency at the ROI, use: `lookup(signal, best_t12ESum)`
effbb = lookup(signal, best_t12ESum)

neutron = restrict(get_bkg_counts_1D(background[end]), 1000, 2900) |> integral
other = expBkgESum - neutron
# Plotting:

# 2D tHalf map:
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
      b  = $(round(expBkgESum, sigdigits = 3)) 
      T12 ≥  $(round(ThalfbbESum, sigdigits=3)) yr 
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

t12(t, e, b) = get_tHalf(
    SNparams["W"],
    SNparams["foilMass"],
    SNparams["Nₐ"],
    t,
    SNparams["a"],
    e,
    (b)/ SNparams["tYear"] * t,
    α;
    approximate="table"
)

let 
    f = Figure(size=(600, 400))
    a = Axis(f[1,1], xlabel = "running time (yr)", ylabel = "sensitivity (yr)", limits= (0,6, nothing, nothing))
    p = lines!(a, t, t12.(t, effbb,expBkgESum), label = "FC: $(analysisDict[:Eres]) resolution; B: off, TKReconstruct \nneutron config = $(background[end].histName)")
    axislegend(a, position = :lt)
    saveName = savename("sensitivity_in_time", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f
end


# Total background model:
bkg_hists = get_bkg_counts_1D.(background)

with_theme(theme_latexfonts()) do
    f = Figure()
    ax = Axis(
        f[1,1], 
        xlabel = analysisDict[:mode], 
        ylabel = "counts", 
        yscale = log10, 
        limits = (300, 3500, 1e-5, 1e6),
        title = "Total background model"
    )
    
    colors = colorschemes[:tol_bright]
    labels = [b.isotopeName for b in background]
    labels[end] = "neutron_external\n$(analysisDict[:neutron_config])"
    st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
	for i=2:length(bkg_hists)
		hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
	end
    
    ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
    ax.xticks = 0:500:3500
    Legend(f[2,1], ax, orientation=:horizontal, fontsize=8, nbanks = 3)
    saveName = savename("background_model", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f
end

# Background by source:
ROI_a, ROI_b = best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]

bkgs = [sum(bincounts(restrict(b, ROI_a, ROI_b)))  for b in get_bkg_counts_1D.(background)]

let 
    saveName = savename("background_counts_ROI=$(ROI_a):$(ROI_b)_", analysisDict, "md")
    open(plotsdir("backgroundTables", saveName), "w") do f
        labels = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_radon", "Tl208_foil_bulk", "neutron_external\n$(analysisDict[:neutron_config])", "total"]
        header = ["process", "bkg counts in ROI"]
        pretty_table(f,
            DataFrame(
                # process = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_radon", "Tl208_foil_bulk", "K40_foil_bulk", "Pa234m_foil_bulk", "neutron_external\n$(analysisDict[:neutron_config])", "total"],
                process = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_radon", "Tl208_foil_bulk", "neutron_external_$(analysisDict[:neutron_config])", "total"],
                counts = vcat(bkgs, sum(bkgs)),
            ),
            header = header,
            alignment = [:l, :l],
            backend = Val(:markdown),
        )
   end
end


# save signal tables
function save_sensitivity_table(
    signals, 
    background, 
    outDir;
    analysisDict = analysisDict,
    α = α,
    SNparams = SNparams,
    )
    df = DataFrame(signal = [], ROI = [], bkg_counts = [], eff = [], t12 = [])
    for s in signals
        t12MapESum = get_tHalf_map(SNparams, α, s, background...; approximate ="formula")
        best_t12ESum = get_max_bin(t12MapESum)
        expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
        effbb = lookup(s, best_t12ESum)
        ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)
        push!(df, (s.isotopeName, "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV", round(expBkgESum, sigdigits = 3), round(effbb, sigdigits = 3), round(ThalfbbESum, sigdigits=3)))
    end

    saveName = savename("sensitivity_", analysisDict, "md")
    open(plotsdir(outDir, saveName), "w") do f
        header = ["signal", "ROI", "bkg counts in ROI", "efficiency", "t12"]
        pretty_table(f, df, header = header, backend = Val(:markdown))
    end
end

signals = [
    get_process("bb0nu_foil_bulk", data_processes),
    get_process("bb0nuM1_foil_bulk", data_processes),
    get_process("bb0nuM2_foil_bulk", data_processes)
]

set_nTotalSim!.( signals, 1e8 )

save_sensitivity_table(signals, background, "sensitivityTables")


#save background tables for each signals ROI
function save_background_table(
    signal, 
    background, 
    analysisDict = analysisDict,
    α = α,
    SNparams = SNparams,
    )
    t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
    best_t12ESum = get_max_bin(t12MapESum)


    ROI_a, ROI_b = best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]

    bkgs = [sum(bincounts(restrict(b, ROI_a, ROI_b)))  for b in get_bkg_counts_1D.(background)]

    analysisDict[:signal] = signal.isotopeName

    saveName = savename("background_counts", analysisDict, "md")
    open(plotsdir("backgroundTables", saveName), "w") do f
        process = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_radon", "Tl208_foil_bulk", "K40_foil_bulk", "Pa234m_foil_bulk", "neutron_external\n$(analysisDict[:neutron_config])", "total"]
        header = ["process", "bkg counts in ROI"]
        pretty_table(f,
            DataFrame(
                process = process,
                counts = vcat(bkgs, sum(bkgs)),
            ),
            header = header,
            alignment = [:l, :l],
            backend = Val(:markdown),
        )
    end
end

for s in signals
    save_background_table(s, background, "backgroundTables")
end



function get_sensitivities_vs_time(
        signal,
        background,
        SNparams;
        neutron_bkg = 0.0,
        effFactor = 1.0
    )
    t = range(0, 5, 100)
    sensitivities = []
    
    t12(t, e, b) = get_tHalf(
        SNparams["W"],
        SNparams["foilMass"],
        SNparams["Nₐ"],
        t,
        SNparams["a"],
        e*effFactor,
        (b+neutron_bkg)/ SNparams["tYear"] * t,
        α;
        approximate="table"
    )
    t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
    best_t12ESum = get_max_bin(t12MapESum)
    expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
    effbb = lookup(signal, best_t12ESum)
    append!(sensitivities, t12.(t, effbb,expBkgESum))
    return sensitivities
end


#### bb0nuM1_foil_bul
signal = get_process("bb0nuM1_foil_bulk", data_processes)
set_nTotalSim!( signal, 1e8 )

background[end] =  get_process("neutron_external", hist_processes, "current shielding")

t_current_nu0M1 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams
)



#################

background[end] = get_process("neutron_external", hist_processes, "full shielding")

t_full_nu0M1 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams
)



let
    f = Figure(size=(600, 400))
    a = Axis(
        f[1,1], 
        xlabel = "running time (yr)", 
        ylabel = "sensitivity (yr)", 
        limits= (0,5, nothing, nothing),
        title = L"Sensitivity for $0\nu\beta\beta\chi^0$ at 90% CL",
        # yscale = log10
        )
    p = lines!(a, t, t_current_nu0M1, label = "current shielding", linewidth = 4, linestyle=(:dash, :loose))
    lines!(a, t, t_full_nu0M1, label = "full shielding no floor flux", linewidth = 4, linestyle=(:dot, :loose))
    hlines!(a, [1.2e23], color = :black, linestyle = :dash, label = "CUPID-0 90% CL", linewidth = 2)
    axislegend(a, position = :lt, patchsize = (30, 20))
    saveName = savename("sensitivity_in_time_nu0M1", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f 
end


#### M2
signal = get_process("bb0nuM2_foil_bulk", data_processes)
background[end] = get_process("neutron_external", hist_processes, "current shielding")

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!(signal, 1e8 )

t_current_nu0M2 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams
)



#################
background[end] = get_process("neutron_external", hist_processes, "full shielding")

t_full_nu0M2 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams
)


let
    f = Figure(size=(600, 400))
    a = Axis(
        f[1,1], 
        xlabel = "running time (yr)", 
        ylabel = "sensitivity (yr)", 
        limits= (0,5, nothing, nothing),
        title = L"Sensitivity for $0\nu\beta\beta\chi^0\chi^0$ at 90% CL"
        )
    p = lines!(a, t, t_current_nu0M2, label = "current shielding", linewidth = 4, linestyle=(:dash, :loose))
    p = lines!(a, t, t_full_nu0M2, label = "full shielding no floor flux", linewidth = 4, linestyle=(:dot, :loose))
    hlines!(a, [1.4e22], color = :black, linestyle = :dash, label = "CUPID-0 90% CL", linewidth = 2)
    axislegend(a, position = :lt, patchsize = (30, 20))
    saveName = savename("sensitivity_in_time_nu0M2", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f 
end


#### 0nu
signal = get_process("bb0nu_foil_bulk", data_processes)
set_nTotalSim!( signal, 0.98e8 )


background[end] = get_process("neutron_external", hist_processes, "current shielding")

t_current_nu0 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams
)


#################
background[end] = get_process("neutron_external", hist_processes, "full shielding")

t_full_nu0 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams
)

t_none100_nu0 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    neutron_bkg = 0.01*100
)


t_none500_nu0 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    neutron_bkg = 0.01*500
)

let
    f = Figure(size=(600, 400))
    a = Axis(
        f[1,1], 
        xlabel = "running time (yr)", 
        ylabel = "sensitivity (yr)", 
        limits= (0,5, nothing, nothing),
        title = L"Sensitivity for $0\nu\beta\beta$ at 90% CL"
        )
    p = lines!(a, t, t_current_nu0, label = "current shielding", linewidth = 4)
    lines!(a, t, t_full_nu0, label = "full shielding", linewidth = 4)
    lines!(a, t, t_none100_nu0, label = "no neutron shielding (x100)", linewidth = 4)
    lines!(a, t, t_none500_nu0, label = "no neutron shielding (x500)", linewidth = 4)
    hlines!(a, [4.6e24], color = :black, linestyle = :dash, label = "CUPID-0 90% CL", linewidth = 2)
    axislegend(a, position = :lt, patchsize = (30, 20))
    saveName = savename("sensitivity_in_time_nu0", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f 
end


## RH L
signal = get_process("bb0nu_foil_bulk", data_processes)
set_nTotalSim!( signal, 0.98e8 )

background[end] = get_process("neutron_external", hist_processes, "full shielding")

t_full_RH_L = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    effFactor = 0.489
)

background[end] = get_process("neutron_external", hist_processes, "current shielding")

t_curr_RH_L = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    effFactor = 0.489
)

let
    f = Figure(size=(600, 400))
    a = Axis(
        f[1,1], 
        xlabel = "running time (yr)", 
        ylabel = "sensitivity (yr)", 
        limits= (0,5, nothing, nothing),
        title = L"Sensitivity for $0\nu\beta\beta$ $\lambda$ (V+A) at 90% CL"
        )
    p = lines!(a, t, t_full_RH_L, label = L"full shielding $$", linewidth = 4)
    lines!(a, t, t_curr_RH_L, label = L"current shielding $$", linewidth = 4)
    hlines!(a, [1.6e23], color = :black, linestyle = :dash, label = L"best $^{82}$Se: $\langle \lambda \rangle$", linewidth = 2)
    band!([0,5], [4.58e23], [13.35e23], color = (:red, 0.4), label = L"best world: $\langle \lambda \rangle$")
    
    # hlines!(a, [2.2e23], color = :red, linestyle = :dash, label = L"$\lambda$: NEMO3", linewidth = 2)
    axislegend(a, position = :lt, patchsize = (30, 20))
    saveName = savename("sensitivity_in_time_nu0_V+A_L", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f 
end

## RH E

signal = get_process("bb0nu_foil_bulk", data_processes)
set_nTotalSim!( signal, 0.98e8 )

background[end] = get_process("neutron_external", hist_processes, "full shielding")

t_full_RH_e = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    effFactor = 0.888
)

background[end] = get_process("neutron_external", hist_processes, "current shielding")

t_curr_RH_e = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    effFactor = 0.888
)

let
    f = Figure(size=(600, 400))
    a = Axis(
        f[1,1], 
        xlabel = "running time (yr)", 
        ylabel = "sensitivity (yr)", 
        limits= (0,5, nothing, nothing),
        title = L"Sensitivity for $0\nu\beta\beta$ $\eta$ (V+A) at 90% CL"
        )
    p = lines!(a, t, t_full_RH_e, label = L"full shielding $$", linewidth = 4)
    lines!(a, t, t_curr_RH_e, label = L"current shielding $$", linewidth = 4)
    # hlines!(a, [1.6e23], color = :black, linestyle = :dash, label = L"$\lambda$: NEMO3", linewidth = 2)
    hlines!(a, [2.2e23], color = :black, linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)
    band!([0,5], [9.93e23], [38.81e23], color = (:red, 0.4), label = L"best world: $\langle \eta \rangle$")
    axislegend(a, position = :lt, patchsize = (30, 20))
    saveName = savename("sensitivity_in_time_nu0_V+A_e", analysisDict, "png")
    safesave(plotsdir("example", analysisDict[:mode], saveName), f, px_per_unit = 6)
    f 
end

