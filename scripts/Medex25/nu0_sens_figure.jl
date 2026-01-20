using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, FHist, ColorSchemes


include(scriptsdir("ND/results/best_rois.jl"))
include(srcdir("params/Params.jl"))

analysisDict    = Dict(:signal => "bb0nu_foil_bulk", :bining => (300, 3500), :bin_width => 100, :mode => "sumE")

roi             = bb0nu_roi
roi[:sumE]      = (300, 3500)
Bin_low, Bin_high, bin_width = analysisDict[:bining][1], analysisDict[:bining][2], analysisDict[:bin_width]
# backgrounds     = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_wire_surface", "Tl208_foil_bulk", "K40_foil_bulk", "Pa234m_foil_bulk", "gamma_experimental_surface"]
backgrounds     = ["bb_foil_bulk", "Bi214_foil_bulk", "Bi214_wire_surface", "Tl208_foil_bulk", "gamma_experimental_surface"]

all_processes   = load_data_processes(datadir("mva/fal5_TKrec_J40"), analysisDict[:mode], fwhm=0.0, roi = roi)
signal          = get_process(analysisDict[:signal], all_processes) |> first
background      = [get_process(b, all_processes) |> first for b in backgrounds]

colors          = ["#003865", "#FFB948", "#52473B", "#9A3A06", ]

set_nTotalSim!(signal, 0.1e8)
set_nTotalSim!(background[1], 1e8)
set_activity!(background[end], 0.185)


# set_activity!(background[1], 0.0019235366786346892)
for b in background
    set_bins!(b, Bin_low:bin_width:Bin_high+bin_width)
end
set_bins!(signal, Bin_low:bin_width:Bin_high+bin_width)


set_signal!(background[1], false)

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


nu0_t_sens = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    neutron_bkg = 0.0,
    effFactor = 1.0
)

get_bkg_counts_ROI(best_t12ESum, background...)

sum(get_bkg_counts_ROI(best_t12ESum, background[i]) for i in 1:length(background))
    
let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e23, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{0\nu\beta\beta}$ decay",)


    h1 = hlines!(ax, [4.6e24], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, nu0_t_sens, color = colors[1], linewidth = 5.5)
    lines!(ax2, t |> collect, nu0_t_sens, color = colors[1], linewidth = 5.5)
    s1 = scatter!(ax, [17.5], [4.64e24], marker = :circle, markersize = 22, color = colors[2])
    
    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [p1, s1, h1], 
        [
            L"$$  frequentist limit @ 90% CL", 
            L"$$  bayesian limit @ 90% CI",
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (40, 50),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity",
        )
    save(scriptsdir("Medex25", "nu0_sens_figure.png"), f, px_per_unit = 5)
    f
    
end

signalM1 = get_process("bb0nuM1_foil_bulk", all_processes) |> first
set_nTotalSim!(signalM1, 0.1e8)
set_bins!(signalM1, Bin_low:bin_width:Bin_high+bin_width)

nu0M1_t_sens = get_sensitivities_vs_time(
    signalM1,
    background,
    SNparams
)

signalM2 = get_process("bb0nuM2_foil_bulk", all_processes) |> first
set_nTotalSim!(signalM2, 0.1e8)
set_bins!(signalM2, Bin_low:bin_width:Bin_high+bin_width)
nu0M2_t_sens = get_sensitivities_vs_time(
    signalM2,
    background,
    SNparams
)

signalRH20 = get_process("RH020_foil_bulk", all_processes) |> first
set_nTotalSim!(signalRH20, 0.1e8)
set_bins!(signalRH20, Bin_low:bin_width:Bin_high+bin_width)
nu0RH20_t_sens = get_sensitivities_vs_time(
    signalRH20,
    background,
    SNparams
)

signalRH50 = get_process("RH050_foil_bulk", all_processes) |> first
set_nTotalSim!(signalRH50, 0.1e8)
set_bins!(signalRH50, Bin_low:bin_width:Bin_high+bin_width)
nu0RH50_t_sens = get_sensitivities_vs_time(
    signalRH50,
    background,
    SNparams
)

let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e21, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e21, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to BSM decays $$",)

    hnu0 = hlines!(ax, [4.6e24], color = colors[1], linestyle = (:dash, :dense), linewidth = 3.5)
    hnu0M1 = hlines!(ax, [1.2e23], color = colors[2], linestyle = (:dash, :dense), linewidth = 3.5)
    hnu0M2 = hlines!(ax, [1.4e22], color = colors[3], linestyle = (:dash, :dense), linewidth = 3.5)

    pnu0 = lines!(ax, exp, nu0_t_sens, color = colors[1], linewidth = 5.5)
    lines!(ax2, t |> collect, nu0_t_sens, color = colors[1], linewidth = 5.5)
    pnu0M1 = lines!(ax, exp, nu0M1_t_sens, color = colors[2], linewidth = 5.5)
    pnu0M2 = lines!(ax, exp, nu0M2_t_sens, color = colors[3], linewidth = 5.5)
    pRH020 = lines!(ax, exp, nu0RH20_t_sens, color = (colors[4], 1), linewidth = 5.5)
    pRH050 = lines!(ax, exp, nu0RH50_t_sens, color = (colors[4], 0.5), linewidth = 5.5)


    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    Legend(f[1,2], 
        [pnu0, hnu0, pnu0M1, hnu0M1, pnu0M2, hnu0M2, pRH020, pRH050], 
        [
            L"  $0\nu\beta\beta$", 
            L"  $$CUPID-0 @ 90% CI",
            L"  $0\nu\beta\beta\chi^0$, n=1",
            L"  $$CUPID-0 @ 90% CI",
            L"  $0\nu\beta\beta\chi^0\chi^0$, n=3",
            L"  $$CUPID-0 @ 90% CI",
            L"  $\nu_R\nu_L\beta\beta; ~\Kappa^{2\nu}_\epsilon = 0.20$",
            L"  $\nu_R\nu_L\beta\beta; ~\Kappa^{2\nu}_\epsilon = 0.50$",
        ], 
        patchsize = (40, 42),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity",
    )
    save(scriptsdir("Medex25", "nu0_bsm_sens_figure.png"), f, px_per_unit = 5)
    f
end





