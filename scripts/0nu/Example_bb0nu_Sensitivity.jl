using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes, SensitivityModule, CairoMakie, UnROOT, StatsBase
using LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta
using Revise
Revise.track(SensitivityModule)

#save background tables for each signals ROI
function save_background_table(
    signal, 
    background,
    outDir;
    labels = labels,
    analysisDict = analysisDict,
    α = α,
    SNparams = SNparams,
    ROI = nothing,
)
    t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
    best_t12ESum = get_max_bin(t12MapESum)

    if( ROI == nothing)
        ROI_a, ROI_b = best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]
    else
        ROI_a, ROI_b = ROI
    end

    
    bkgs = [sum(bincounts(restrict(b, ROI_a, ROI_b))) ± sqrt(sum(sumw2(restrict(b, ROI_a, ROI_b)))) for b in get_bkg_counts_1D.(background)]

    analysisDict[:signal] = signal.isotopeName

    saveName = savename("background_counts", analysisDict, "md")
    open(plotsdir(outDir, saveName), "w") do f
        header = ["process", "bkg counts in ROI"]
        pretty_table(f,
            DataFrame(
                process = vcat(labels, "total"),
                counts = vcat(bkgs, sum(bkgs)),
            ),
            header = header,
            alignment = [:l, :l],
            backend = Val(:markdown),
        )
    end
end


# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
include(srcdir("params/Params.jl"))

# Dictionary with the analysis parameters. (Nice for when you want to save stuff and make sure you know what params you used in this analysis.)
analysisDict = Dict(
    :Bfield => "Boff", # magnetic field on/off
    :Eres => "8perc", # FWHM of the OMs (sorry for the naming...)
    :mode => "sumE", 
    :trackAlgo => "TKrec",
    :signal => "bb0nu",
    :neutron_config => "current_shielding"
)

# files_directory = "fal5_$(analysisDict[:Eres])_$(analysisDict[:Bfield])_$(analysisDict[:trackAlgo])_twoDistinct_edep_bcu"
files_directory = "fal5_8perc_Boff_TKrec_evis_bcu_J38"

# Load all the processes in the directory. Function `load_processes` takes two arguments:
#  1. dir::String -> the name of the directory where the root files are stored
#  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
data_processes = load_data_processes(
    files_directory, 
    analysisDict[:mode],
    fwhm = 0.0
)

hist_processes = load_hist_processes(
    files_directory,  
    analysisDict[:mode]
)

for c in ["current_shielding"]
# begin
    analysisDict[:neutron_config] = c
    # declare which process is signal
    global signal = get_process("bb0nu_foil_bulk", data_processes) |> first
    # signal = get_process("bb0nuM1_foil_bulk", data_processes)
    # signal = get_process("bb0nuM2_foil_bulk", data_processes)

    # declare background processes
    global background = [
        get_process("bb_foil_bulk", data_processes) |> first,
        get_process("Bi214_foil_bulk", data_processes) |> first,
        get_process("Bi214_wire_surface", data_processes) |> first,
        get_process("Tl208_foil_bulk", data_processes) |> first,
        # get_process("K40_foil_bulk", data_processes) |> first,
        # get_process("Pa234m_foil_bulk", data_processes) |> first,
        get_process("neutron_external", hist_processes, analysisDict[:neutron_config])
        # get_process("Bi214_PMT_glass_bulk", data_processes) |> first
    ]

    labels = [L"$2\nu\beta\beta$", L"$^{214}$Bi",L"Radon $$", L"$^{208}$Tl", "external"]
    # labels = [L"$2\nu\beta\beta$", L"$^{214}$Bi", L"Radon $$", L"$^{208}$Tl", L"$^{40}$K", L"$^{234m}$Pa", "neutrons (5-sided)"]
    # labels = [L"$2\nu\beta\beta$", L"$^{214}$Bi", L"Radon $$", L"$^{208}$Tl"]

    # if(analysisDict[:neutron_config] == "full_shielding")
    #     labels[end] = "neutrons (6-sided)"
    # elseif (analysisDict[:neutron_config] == "iron_shielding")
    #     labels[end] = "neutrons (0-sided)"
    # elseif (analysisDict[:neutron_config] == "no_french_wall_shielding")
    #     labels[end] = "neutrons (5-sided)"
    # elseif (analysisDict[:neutron_config] == "italian_shielding")
    #     labels[end] = "neutrons (5-sided)"
    # elseif (analysisDict[:neutron_config] == "current_shielding")
    #     labels[end] = "neutrons"
    # else
    #     labels[end] = "neutrons"
    # end

    # set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
    set_signal!(background[1], false)

    # set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
    # set_nTotalSim!( signal, 0.98e8 )
    set_nTotalSim!( signal, 1e8 )
    # set_nTotalSim!( signal, 1e8 )
    # set_nTotalSim!( background[1], 0.99e8 )
    # set_nTotalSim!( background[2], 0.96e8 )
    # set_nTotalSim!( background[3], 1e8 )
    # set_nTotalSim!( background[4], 0.76e8 )
    # set_nTotalSim!( background[5], 1e8 )
    # set_nTotalSim!( background[6], 1e8 )

    set_nTotalSim!( background[1], 1e8 )
    set_nTotalSim!( background[2], 1e8 )
    set_nTotalSim!( background[3], 1e8 )
    set_nTotalSim!( background[4], 1e8 )
    # set_nTotalSim!( background[5], 1e8 )
    # set_nTotalSim!( background[6], 1e8 )

    println("Processes initialized.")

    # For now I only implemented 90% CL estimation, maybe someday I'll get to more...
    global α = 1.64485362695147

    global t = range(0, 5, 100)

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
    @show best_t12ESum = get_max_bin(t12MapESum)

    # If you want additional info, like background counts in the ROI, use: `get_bkg_counts_ROI`.
    @show expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)

    # To get the signal efficiency at the ROI, use: `lookup(signal, best_t12ESum)`
    @show effbb = lookup(signal, best_t12ESum)

    neutron = restrict(get_bkg_counts_1D(background[end]), 2750, 2950) |> integral
    other = expBkgESum - neutron
    # Plotting:

    # 2D tHalf map:
    @show ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

    set_activity!(signal, halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], ThalfbbESum * 365 * 24 * 3600))

    lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
          b  = $(round(expBkgESum, digits = 2)) 
          T12 ≥  $(round(ThalfbbESum, sigdigits=3)) yr 
          ε = $(round(effbb, digits = 2)*100)%"

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


    # # Sensitivity as a function of detector life-time:

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
        # p = lines!(a, t, t12.(t, effbb,expBkgESum), label = "FC: $(analysisDict[:Eres]) resolution; B: off, TKReconstruct \nneutron config = $(background[end].histName)")
        p = lines!(a, t, t12.(t, effbb,expBkgESum) , label = "SuperNEMO sensitivity estimate", linewidth = 4)
        axislegend(a, position = :lt)
        saveName = savename("sensitivity_in_time", analysisDict, "png")
        safesave(plotsdir("for_Will", analysisDict[:mode], saveName), f, px_per_unit = 6)
        f
    end


    # Total background model:
    bkg_hists = get_bkg_counts_1D.(background)
    sig_hist = get_bkg_counts_1D(signal)


    with_theme(theme_latexfonts()) do
        f = Figure()
        ax = Axis(
            f[1,1], 
            # xlabel = analysisDict[:mode], 
            xlabel = L"$E_1 + E_2$ (keV)", 
            ylabel = L"counts / $17.5$kg.yr exposure / $100$ keV", 
            yscale = log10, 
            limits = (300, 3500, 1e-4, 1e6),
            # limits = (2500, 3500, 0, 3),
            title = "Total background model\nsummed 2-electron energy"
        )
        
        colors = colorschemes[:Set2_7]
        # labels = [b.isotopeName for b in background]
        # labels[end] = "neutron_external\n$(analysisDict[:neutron_config])"
        st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
        errorbars!(ax, sum(bkg_hists), color = :black, whiskerwidth = 7)
        
        for i=2:length(bkg_hists)
            hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
            # errorbars!(ax, sum(bkg_hists[i:end]), color = :black, whiskerwidth = 7)
            
        end
        # lines!(ax, midpoints(binedges(sig_hist)), bincounts(sig_hist), label = signal.isotopeName, color = :red, linestyle = :dash, linewidth = 2.5)
        # errorbars!(ax, sum(bkg_hists), color = :black, whiskerwidth = 7)

        ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
        ax.xticks = 0:500:3500
        Legend(f[2,1], ax, orientation=:horizontal, fontsize=8, nbanks = 2)
        saveName = savename("background_model", analysisDict, "png")
        wsave(plotsdir("for_Will", analysisDict[:mode], saveName), f, px_per_unit = 6)
        f
    end


    # # save signal tables
    # function save_sensitivity_table(
    #     signals, 
    #     background, 
    #     outDir;
    #     analysisDict = analysisDict,
    #     α = α,
    #     SNparams = SNparams,
    #     )
    #     df = DataFrame(signal = [], ROI = [], bkg_counts = [], eff = [], t12 = [])
    #     for s in signals
    #         t12MapESum = get_tHalf_map(SNparams, α, s, background...; approximate ="table")
    #         best_t12ESum = get_max_bin(t12MapESum)
    #         bkg_hists = get_bkg_counts_1D.(background)
    #         b = integral(restrict(sum(bkg_hists), best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]))
    #         Δb = sqrt(sum(sumw2(restrict(sum(bkg_hists), best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]))))
    #         expBkgESum = b ± Δb
    #         effbb = lookup(s, best_t12ESum)
    #         ThalfbbESum = round(get_tHalf(SNparams, effbb, Measurements.value(expBkgESum), α), digits=2)
    #         push!(df, (s.isotopeName, "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV", round(expBkgESum, digits = 2), round(effbb, digits = 2), round(ThalfbbESum, sigdigits=3)))
    #     end

    #     saveName = savename("sensitivity_", analysisDict, "md")
    #     open(plotsdir(outDir, saveName), "w") do f
    #         header = ["signal", "ROI", "bkg counts in ROI", "efficiency", "t12"]
    #         pretty_table(f, df, header = header, backend = Val(:markdown))
    #     end
    # end

    signals = [
        get_process("bb0nu_foil_bulk", data_processes) |> first,
        # get_process("bb0nuM1_foil_bulk", data_processes),
        # get_process("bb0nuM2_foil_bulk", data_processes)
    ]

    # set_nTotalSim!.( signals[2:end], 1e8 )

    # save_sensitivity_table(signals, background, "LSM_report/sensitivityTables/sumE")

    for s in signals
        # save_background_table(s, background, "for_Will/backgroundTables/sumE"; labels=labels, ROI = (0, 3500))
        save_background_table(s, background, "for_Will/backgroundTables/sumE"; labels=labels)
    end
end;

# begin
    global α = 1.64485362695147
    global t = range(0, 5, 100)

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


    ### bb0nuM1_foil_bul
    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface])

    signal = get_process("bb0nuM1_foil_bulk", data_processes) |> first
    set_nTotalSim!( signal, 1e8 )

    # background[end] =  get_process("neutron_external", hist_processes, "current_shielding")

    t_current_nu0M1 = get_sensitivities_vs_time(
        signal,
        background,
        SNparams
    )

    # #################
    # background[end] = get_process("neutron_external", hist_processes, "iron_shielding")

    # t_iron_nu0M1 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # #################
    # background[end] = get_process("neutron_external", hist_processes, "no_french_wall_shielding")

    # t_italian_nu0M1 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # #################

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_nu0M1 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface]*100)

    # signal = get_process("bb0nuM1_foil_bulk", data_processes)
    # set_nTotalSim!( signal, 1e8 )

    # background[end] =  get_process("neutron_external", hist_processes, "current_shielding")

    # t_current_nu0M1_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_nu0M1_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # let
    #     f = Figure(size=(600, 400))
    #     a = Axis(
    #         f[1,1], 
    #         xlabel = "running time (yr)", 
    #         ylabel = "sensitivity (yr)", 
    #         limits= (0,5, nothing, nothing),
    #         title = L"Sensitivity for $0\nu\beta\beta\chi^0$ at 90% CL",
    #         # yscale = log10
    #         )
    #     p1 = lines!(a, t, t_full_nu0M1, label = "6-sided", linewidth = 2.5)#, linestyle=(:dot, :dense))
    #     # lines!(a, t, t_italian_nu0M1, label = "5-sided", linewidth = 2.5)#, linestyle=(:dashdotdot, :dense))
    #     p2 = lines!(a, t, t_current_nu0M1, label = "4-sided", linewidth = 2.5)#, linestyle=(:dash, :dense))
    #     # lines!(a, t, t_iron_nu0M1, label = "0-sided", linewidth = 2.5)#, linestyle=(:dashdot, :dense))
    #     p3 = lines!(a, t, t_full_nu0M1_wo, label = "6-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[1])
    #     p4 =lines!(a, t, t_current_nu0M1_wo, label = "4-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[2])
    #     h= hlines!(a, [1.2e23], color = :black, linestyle = :solid, label = "CUPID-0 90% CL", linewidth = 2)
    #     axislegend(a, position = :rb, patchsize = (20, 7), patchlabelgap = 6, nbanks = 2)
    #         # a.yticks = ([1e22, 1e23, 1e24], [L"10^{22}", L"10^{23}", L"10^{24}"])
    #     saveName = savename("sensitivity_in_time_nu0M1", analysisDict, "png")
    #     safesave(plotsdir("LSM_report", "sensitivity_over_time", analysisDict[:mode], saveName), f, px_per_unit = 6)
    #     f 
    # end


    # #### M2
    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface])
    signal = get_process("bb0nuM2_foil_bulk", data_processes) |> first
    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    # # set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
    set_nTotalSim!(signal, 1e8 )

    t_current_nu0M2 = get_sensitivities_vs_time(
        signal,
        background,
        SNparams
    )

    # #################

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_nu0M2 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface]*100.0)

    # signal = get_process("bb0nuM2_foil_bulk", data_processes)
    # set_nTotalSim!( signal, 1e8 )

    # background[end] =  get_process("neutron_external", hist_processes, "current_shielding")

    # t_current_nu0M2_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_nu0M2_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # #################
    # background[end] = get_process("neutron_external", hist_processes, "no_french_wall_shielding")

    # t_italian_nu0M2 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )


    # #################
    # background[end] = get_process("neutron_external", hist_processes, "iron_shielding")

    # t_iron_nu0M2 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # let
    #     f = Figure(size=(600, 400))
    #     a = Axis(
    #         f[1,1], 
    #         xlabel = "running time (yr)", 
    #         ylabel = "sensitivity (yr)", 
    #         limits= (0,5, nothing,nothing),
    #         # yscale= log10,
    #         title = L"Sensitivity for $0\nu\beta\beta\chi^0\chi^0$ at 90% CL"
    #         )
    #     p1 = lines!(a, t, t_full_nu0M2, label = "6-sided", linewidth = 2.5)#, linestyle=(:dot, :dense))
    #     # lines!(a, t, t_italian_nu0M2, label = "5-sided", linewidth = 2.5)#, linestyle=(:dashdotdot, :dense))
    #     p2 = lines!(a, t, t_current_nu0M2, label = "4-sided", linewidth = 2.5)#, linestyle=(:dash, :dense))
    #     # lines!(a, t, t_iron_nu0M2, label = "0-sided", linewidth = 2.5)#, linestyle=(:dashdot, :dense))
    #     p3 = lines!(a, t, t_full_nu0M2_wo, label = "6-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[1])
    #     p4 =lines!(a, t, t_current_nu0M2_wo, label = "4-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[2])
    #     hlines!(a, [1.4e22], color = :black, linestyle = :solid, label = "CUPID-0 90% CL", linewidth = 2)
    #     axislegend(a, position = :rb, patchsize = (20, 7), patchlabelgap = 6, nbanks = 2)
    #     # a.yticks = ([1e22, 1e23, 1e24], [L"10^{22}", L"10^{23}", L"10^{24}"])
    #     saveName = savename("sensitivity_in_time_nu0M2", analysisDict, "png")
    #     safesave(plotsdir("LSM_report","sensitivity_over_time", analysisDict[:mode], saveName), f, px_per_unit = 6)
    #     f 
    # end


    #### 0nu
    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface])

    signal = get_process("bb0nu_foil_bulk", data_processes) |> first
    set_nTotalSim!( signal, 1e8 )
    # set_nTotalSim!( signal, 0.98e8 )


    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    t_current_nu0 = get_sensitivities_vs_time(
        signal,
        background,
        SNparams
    )

    # background[end] = get_process("neutron_external", hist_processes, "iron_shielding")

    # t_iron_nu0 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # background[end] = get_process("neutron_external", hist_processes, "no_french_wall_shielding")

    # t_italian_nu0 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # #################
    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_nu0 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    ### Adjust radon activity
    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface]*100.0)
    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    # t_current_nu0_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_nu0_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams
    # )

    # t_none100_nu0 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     neutron_bkg = 0.01*100
    # )


    # t_none500_nu0 = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     neutron_bkg = 0.01*500
    # )

    # let
    #     f = Figure(size=(700, 500), fontsize = 16)
    #     a = Axis(
    #         f[1,1], 
    #         xlabel = "exposure (kg.yr)", 
    #         ylabel = "sensitivity (yr)", 
    #         limits= (0,5*6.11, nothing, nothing),
    #         title = L"Sensitivity for $0\nu\beta\beta$ at 90% CL"
    #         )
    #     # p1 = lines!(a, t, t_full_nu0, label = "6-sided", linewidth = 2.5)#, linestyle=(:dot, :dense))
    #     # lines!(a, t, t_italian_nu0, label = "5-sided", linewidth = 2.5)#, linestyle=(:dashdotdot, :dense))
    #     # p2 = lines!(a, t, t_current_nu0, label = "4-sided", linewidth = 2.5)#, linestyle=(:dash, :dense))
    #     p2 = lines!(a, t .* 6.11, t_current_nu0, label = "SuperNEMO sensitivity estimate", linewidth = 4.5)#, linestyle=(:dash, :dense))
    #     # lines!(a, t, t_iron_nu0, label = "0-sided", linewidth = 2.5)#, linestyle=(:dashdot, :dense))
    #     # p3 = lines!(a, t, t_full_nu0_wo, label = "6-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[1])
    #     # p4 =lines!(a, t, t_current_nu0_wo, label = "4-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[2])
    #     # h= hlines!(a, [4.6e24], color = :black, linestyle = :solid, label = "CUPID-0 90% CL", linewidth = 2)
    #     axislegend(a, position = :lt, patchsize = (20, 7), patchlabelgap = 6, nbanks = 2)
    #     # Legend(
    #     #     f[1,1], [h], ["CUPID-0 90% CL"], fontsize=12, 
    #     #     orientation=:horizontal, tellheight = false, tellwidth = false,
    #     #     margin = (-150, 0, -350, 10),
    #     # )
    #     # Legend(
    #     #     f[1,1], [p1, p2], ["6-sided", "4-sided"], "With Anti-Rn", fontsize=12, 
    #     #     tellheight = false, tellwidth = false,
    #     #     margin = (150, 0, -300, 10),
    #     # )
    #     # Legend(
    #     #     f[1,1], [p3, p4], ["6-sided", "4-sided"], "Without Anti-Rn", fontsize=12, 
    #     #     tellheight = false, tellwidth = false,
    #     #     margin = (420, 0, -300, 10),
    #     # )
    #     saveName = savename("sensitivity_in_time_nu0", analysisDict, "png")
    #     safesave(plotsdir("for_Will", "sensitivity_over_time", analysisDict[:mode], saveName), f, px_per_unit = 6)
    #     f 
    # end


    ## RH L
    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface])

    signal = get_process("bb0nu_foil_bulk", data_processes) |> first
    set_nTotalSim!( signal, 1e8 )

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_RH_L = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     effFactor = 0.489
    # )

    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    t_current_RH_L = get_sensitivities_vs_time(
        signal,
        background,
        SNparams;
        effFactor = 0.489
    )

    # ### Adjust radon activity
    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface]*100.0)
    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    # t_current_RH_L_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     effFactor = 0.489

    # )

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_RH_L_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     effFactor = 0.489

    # )

    # # background[end] = get_process("neutron_external", hist_processes, "iron_shielding")

    # # t_iron_RH_L = get_sensitivities_vs_time(
    # #     signal,
    # #     background,
    # #     SNparams;
    # #     effFactor = 0.489
    # # )

    # # background[end] = get_process("neutron_external", hist_processes, "no_french_wall_shielding")

    # # t_italian_RH_L = get_sensitivities_vs_time(
    # #     signal,
    # #     background,
    # #     SNparams;
    # #     effFactor = 0.489
    # # )

    # let
    #     f = Figure(size=(600, 400))
    #     a = Axis(
    #         f[1,1], 
    #         xlabel = "running time (yr)", 
    #         ylabel = "sensitivity (yr)", 
    #         limits= (0,5, nothing, 4e24),
    #         title = L"Sensitivity for $0\nu\beta\beta$ $\lambda$ (V+A) at 90% CL"
    #         )
    #     p1 = lines!(a, t, t_full_RH_L, label = "6-sided", linewidth = 2.5)#, linestyle=(:dot, :dense))
    #     # lines!(a, t, t_italian_RH_L, label = "5-sided", linewidth = 2.5)#, linestyle=(:dashdotdot, :dense))
    #     p2 = lines!(a, t, t_current_RH_L, label = "4-sided", linewidth = 2.5)#, linestyle=(:dash, :dense))
    #     # lines!(a, t, t_iron_RH_L, label = "0-sided", linewidth = 2.5)#, linestyle=(:dashdot, :dense))
    #     p3 = lines!(a, t, t_full_RH_L_wo, label = "6-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[1])
    #     p4 =lines!(a, t, t_current_RH_L_wo, label = "4-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[2])
    #     hlines!(a, [1.6e23], color = :black, linestyle = :solid, label = L"best $^{82}$Se: $\langle \lambda \rangle$", linewidth = 2)
    #     band!([0,5], [4.58e23], [13.35e23], color = (:red, 0.4), label = L"best world: $\langle \lambda \rangle$")
        
    #     # hlines!(a, [2.2e23], color = :red, linestyle = :dash, label = L"$\lambda$: NEMO3", linewidth = 2)
    #     axislegend(a, position = :lt, patchsize = (20, 7), patchlabelgap = 6, nbanks = 2)
    #     saveName = savename("sensitivity_in_time_nu0_V+A_L", analysisDict, "png")
    #     safesave(plotsdir("LSM_report", "sensitivity_over_time", analysisDict[:mode], saveName), f, px_per_unit = 6)
    #     f 
    # end

    # # ## RH E

    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface])
    signal = get_process("bb0nu_foil_bulk", data_processes) |> first
    set_nTotalSim!( signal, 1e8 )

    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_RH_e = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     effFactor = 0.888
    # )

    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    t_current_RH_e = get_sensitivities_vs_time(
        signal,
        background,
        SNparams;
        effFactor = 0.888
    )

    # set_activity!(background[3], BkgActivityParams[:Bi214_wire_surface]*100.0)
    # background[end] = get_process("neutron_external", hist_processes, "full_shielding")

    # t_full_RH_e_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     effFactor = 0.888
    # )

    # background[end] = get_process("neutron_external", hist_processes, "current_shielding")

    # t_current_RH_e_wo = get_sensitivities_vs_time(
    #     signal,
    #     background,
    #     SNparams;
    #     effFactor = 0.888
    # )

    # # background[end] = get_process("neutron_external", hist_processes, "iron_shielding")

    # # t_iron_RH_e = get_sensitivities_vs_time(
    # #     signal,
    # #     background,
    # #     SNparams;
    # #     effFactor = 0.888
    # # )

    # # background[end] = get_process("neutron_external", hist_processes, "no_french_wall_shielding")

    # # t_italian_RH_e = get_sensitivities_vs_time(
    # #     signal,
    # #     background,
    # #     SNparams;
    # #     effFactor = 0.888
    # # )

    # let
    #     f = Figure(size=(600, 400))
    #     a = Axis(
    #         f[1,1], 
    #         xlabel = "running time (yr)", 
    #         ylabel = "sensitivity (yr)", 
    #         limits= (0,5, nothing, nothing),
    #         title = L"Sensitivity for $0\nu\beta\beta$ $\eta$ (V+A) at 90% CL"
    #         )
    #     p1 = lines!(a, t, t_full_RH_e, label = "6-sided", linewidth = 2.5)#, linestyle=(:dot, :dense))
    #     # lines!(a, t, t_italian_RH_e, label = "5-sided", linewidth = 2.5)#, linestyle=(:dashdotdot, :dense))
    #     p2 = lines!(a, t, t_current_RH_e, label = "4-sided", linewidth = 2.5)#, linestyle=(:dash, :dense))
    #     # lines!(a, t, t_iron_RH_e, label = "0-sided", linewidth = 2.5)#, linestyle=(:dashdot, :dense))
    #     p3 = lines!(a, t, t_full_RH_e_wo, label = "6-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[1])
    #     p4 =lines!(a, t, t_current_RH_e_wo, label = "4-sided w/o Anti-Rn factory", linewidth = 2.5, linestyle=(:dash, :dense), color = Makie.wong_colors()[2])
    #     hlines!(a, [2.2e23], color = :black, linestyle = :solid, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)
    #     band!([0,5], [9.93e23], [38.81e23], color = (:red, 0.4), label = L"best world: $\langle \eta \rangle$")
    #     axislegend(a, position = :lt, patchsize = (25, 10), patchlabelgap = 10)
    #     saveName = savename("sensitivity_in_time_nu0_V+A_e", analysisDict, "png")
    #     safesave(plotsdir("LSM_report", "sensitivity_over_time", analysisDict[:mode], saveName), f, px_per_unit = 6)
    #     f 
    # end
# end


# analysisDictSingle = Dict(
#     :Bfield => "Boff", # magnetic field on/off
#     :Eres => "8perc", # FWHM of the OMs (sorry for the naming...)
#     :mode => "singleE", 
#     :trackAlgo => "TIT",
#     :signal => "bb0nu",
#     :neutron_config => "no_french_wall_shielding"
# )

# data_processes = load_data_processes(
#     files_directory, 
#     analysisDictSingle[:mode]
# )


# ######### SINGLE ELECTRON SPECTRUM
# for c in ["full_shielding","iron_shielding","no_french_wall_shielding","current_shielding"]
#     analysisDictSingle[:neutron_config] = c

#     labels = [L"$2\nu\beta\beta$", L"$^{214}$Bi", L"Radon $$", L"$^{208}$Tl", L"$^{40}$K", L"$^{234m}$Pa", "neutrons__"]

#     if(analysisDictSingle[:neutron_config] == "full_shielding")
#         labels[end] = "neutrons (6-sided)"
#     elseif (analysisDictSingle[:neutron_config] == "iron_shielding")
#         labels[end] = "neutrons (0-sided)"
#     elseif (analysisDictSingle[:neutron_config] == "no_french_wall_shielding")
#         labels[end] = "neutrons (5-sided)"
#     elseif (analysisDictSingle[:neutron_config] == "italian_shielding")
#         labels[end] = "neutrons (5-sided)"
#     elseif (analysisDictSingle[:neutron_config] == "current_shielding")
#         labels[end] = "neutrons (4-sided)"
#     else
#         labels[end] = "neutrons"
#     end

#     files_directory = "fal5_$(analysisDictSingle[:Eres])_$(analysisDictSingle[:Bfield])_$(analysisDictSingle[:trackAlgo])_twoDistinct_edep_bcu"

#     # Load all the processes in the directory. Function `load_processes` takes two arguments:
#     #  1. dir::String -> the name of the directory where the root files are stored
#     #  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"

#     fN = ROOTFile("data/sims/neutrons_/neutron_external_single.root")
#     h1 = UnROOT.parseTH(fN[analysisDictSingle[:neutron_config]], raw= false)    

#     hist_processesSingle = HistProcess(
#         h1,
#         analysisDictSingle[:neutron_config],
#         sumEParams[:neutron_external]
#     )

#     # declare which process is signal
#     signal = get_process("bb0nu_foil_bulk", data_processes)
#     # signal = get_process("bb0nuM1_foil_bulk", data_processes)
#     # signal = get_process("bb0nuM2_foil_bulk", data_processes)

#     # declare background processes
#     backgroundSingle = [
#         get_process("bb_foil_bulk", data_processes),
#         get_process("Bi214_foil_bulk", data_processes),
#         get_process("Bi214_wire_surface", data_processes),
#         get_process("Tl208_foil_bulk", data_processes),
#         get_process("K40_foil_bulk", data_processes),
#         get_process("Pa234m_foil_bulk", data_processes),
#         hist_processesSingle
#         # get_process("neutron_external", hist_processesSingle, analysisDictSingle[:neutron_config])
#     ]


#     # set 2nubb to backgroundSingle process (initially it's signal for exotic 2nubb analyses)
#     set_signal!(backgroundSingle[1], false)

#     # set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
#     set_nTotalSim!( signal, 0.98e8 )
#     # set_nTotalSim!( signal, 1e8 )
#     set_nTotalSim!( backgroundSingle[1], 0.99e8 )
#     set_nTotalSim!( backgroundSingle[2], 0.96e8 )
#     set_nTotalSim!( backgroundSingle[3], 1e8 )
#     set_nTotalSim!( backgroundSingle[4], 0.76e8 )
#     set_nTotalSim!( backgroundSingle[5], 1e8 )
#     set_nTotalSim!( backgroundSingle[6], 1e8 )

#     # Total backgroundSingle model:
#     bkg_hists = get_bkg_counts_1D.(backgroundSingle)
#     sig_hist = get_bkg_counts_1D(signal)


#     with_theme(theme_latexfonts()) do
#         f = Figure()
#         ax = Axis(
#             f[1,1], 
#             # xlabel = analysisDictSingle[:mode], 
#             xlabel = L"$E_i$ (keV)", 
#             ylabel = L"counts / $17.5$kg.yr exposure / $100$ keV", 
#             yscale = log10, 
#             limits = (0, 3500, 1e-5, 1e6),
#             # limits = (100, 3500, 0, 3.35e4),
#             title = "Total background model\nsingle-electron energy"
#         )
        
#         colors = colorschemes[:tol_bright]
#         # labels = [b.isotopeName for b in backgroundSingle]
#         # labels[end] = "neutron_external\n$(analysisDictSingle[:neutron_config])"
#         st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black,)
#         for i=2:length(bkg_hists)
#             hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
#         end
#         # lines!(ax, midpoints(binedges(sig_hist)), bincounts(sig_hist), label = signal.isotopeName, color = :red, linestyle = :dash, linewidth = 2.5)
#         errorbars!(ax,  sum(bkg_hists), color = :black, whiskerwidth = 7)
        
#         ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
#         ax.xticks = 0:500:3500
#         Legend(f[2,1], ax, orientation=:horizontal, fontsize=8, nbanks = 2)
#         saveName = savename("background_model", analysisDictSingle, "png")
#         safesave(plotsdir("LSM_report", "background_model", analysisDictSingle[:mode], saveName), f, px_per_unit = 6)
#         f
#     end

#     save_background_table(signal, backgroundSingle, "LSM_report/backgroundTables/singleE/0_1000keV_ROI"; analysisDict = analysisDictSingle, ROI = (0, 1000))

# end



using CSV

dff = DataFrame(time = t, 
                t12_nu0M1 = t_current_nu0M1, 
                t12_nu0M2 = t_current_nu0M2, 
                t12_nu0 = t_current_nu0,
                t12_RH_e = t_current_RH_e,
                t12_RH_L = t_current_RH_L
                )

open("sensitivities.csv", "w") do io
    CSV.write(io, dff)
end


ROI_a, ROI_b = 2700, 3000

bin_counts = [bincounts(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)]
bin_errs = [sqrt.(sumw2(restrict(b, ROI_a, ROI_b))) for b in get_bkg_counts_1D.(background)]
bin_edges = binedges(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b))[1:end-1] |> collect

dfff = DataFrame(
    bin_edges = bin_edges,
    bin_counts_bb = bin_counts[1],
    bin_counts_Bi214 = bin_counts[2],
    bin_counts_Bi214_wire = bin_counts[3],
    bin_counts_Tl208 = bin_counts[4],
    bin_counts_K40 = bin_counts[5],
    bin_counts_Pa234m = bin_counts[6],
    bin_counts_neutrons = bin_counts[7],
    bin_errs_bb = bin_errs[1],
    bin_errs_Bi214 = bin_errs[2],
    bin_errs_Bi214_wire = bin_errs[3],
    bin_errs_Tl208 = bin_errs[4],
    bin_errs_K40 = bin_errs[5],
    bin_errs_Pa234m = bin_errs[6],
    bin_errs_neutrons = bin_errs[7],
)

open("background_counts_per_bin_per_process.csv", "w") do io
    CSV.write(io, dfff)
end

bkgs = [sum(bincounts(restrict(b, ROI_a, ROI_b))) ± sqrt(sum(sumw2(restrict(b, ROI_a, ROI_b)))) for b in get_bkg_counts_1D.(background)]


