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

colors          = ["#003865", "#FFB948", "#52473B", "#9A3A06", "#951272", "#006630", "#005C8A", "#FFB948"]

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
    @show best_t12ESum = get_max_bin(t12MapESum)
    @show expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
    @show effbb = lookup(signal, best_t12ESum)
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

set_activity!(background[3], 20/1000)
nu0_t_sens_Rn20 = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    neutron_bkg = 0.0,
    effFactor = 1.0
)

set_activity!(background[3], 150/1000000)
set_activity!(background[end], 0.185/10)
nu0_t_sens_neutron6_sided = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    neutron_bkg = 0.0,
    effFactor = 1.0
)

# get_bkg_counts_ROI(best_t12ESum, background...)

# sum(get_bkg_counts_ROI(best_t12ESum, background[i]) for i in 1:length(background))
    
let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e23, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{0\nu\beta\beta}$ decay",)


    h1 = hlines!(ax, [4.6e24], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3], 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, nu0_t_sens, color = colors[1], linewidth = 5.5)
    lines!(ax2, t |> collect, nu0_t_sens, color = colors[1], linewidth = 5.5)
    
    s1 = scatter!(ax, [17.5], [4.64e24], marker = :circle, markersize = 22, color = colors[2])
    
    p2 = lines!(ax, exp, nu0_t_sens_Rn20, color = colors[4], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, nu0_t_sens_Rn20, color = colors[4], linewidth = 5.5, linestyle = :dash)

    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [s1, p1, p2, h1], 
        [
            L"$$  bayesian limit @ 90% CI",
            L"$$  frequentist limit @ 90% CL \\Rn @ 0.15~$\textrm{mBq/m^3}$", 
            L"$$  frequentist limit @ 90% CL \\Rn @ 20~$\textrm{mBq/m^3}$", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (40, 70),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure.png"), f, px_per_unit = 5)
    f
    
end


let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e23, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{0\nu\beta\beta}$ decay",)


    h1 = hlines!(ax, [4.6e24], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3], 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, nu0_t_sens, color = colors[1], linewidth = 5.5)
    lines!(ax2, t |> collect, nu0_t_sens, color = colors[1], linewidth = 5.5)
    
    s1 = scatter!(ax, [17.5], [4.64e24], marker = :circle, markersize = 22, color = colors[2])
    
    p2 = lines!(ax, exp, nu0_t_sens_neutron6_sided, color = colors[4], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, nu0_t_sens_neutron6_sided, color = colors[4], linewidth = 5.5, linestyle = :dash)

    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [s1, p1, p2, h1], 
        [
            L"$$  bayesian limit @ 90% CI",
            L"$$  frequentist limit @ 90% CL \\neutron 4-sided", 
            L"$$  frequentist limit @ 90% CL \\neutron 6-sided", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (40, 70),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure_neutron.png"), f, px_per_unit = 5)
    f
    
end


let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e23, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{0\nu\beta\beta}$ decay",)


    h1 = hlines!(ax, [4.6e24], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3], 0.7), linewidth = 1.5)
    
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
        [s1, p1, h1], 
        [
            L"$$  bayesian limit @ 90% CI",
            L"$$  frequentist limit @ 90% CL", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (40, 70),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure_rn150uBq.png"), f, px_per_unit = 5)
    f
    
end


### TH lambda
### For lambda: [9,56-13,35] 1E23 years

t_current_RH_L = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    effFactor = 0.489
)

let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated sensitivity (yr)$$",yscale = log10, limits = (0, maximum(exp), 5e22, 5e24), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{$0\nu\beta\beta \langle \lambda \rangle$ }_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 5e22, 5e24))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta \langle \lambda \rangle$ }$ decay",)



    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3],0.7), linewidth = 1.5)

    p1 = lines!(ax, exp, t_current_RH_L, color = colors[6], linewidth = 5.5)
    lines!(ax2, t |> collect, t_current_RH_L, color = colors[6], linewidth = 5.5)
    
    h1 = hlines!(ax, [1.6e23], color = colors[6], linestyle = :dash, label = L"best $^{82}$Se: $\langle \lambda \rangle$", linewidth = 2)
    b1 = band!(ax, [0,35], [4.58e23], [13.35e23], color = (colors[6], 0.4), label = L"best world: $\langle \lambda \rangle$")
        

    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, h1, b1], 
        [
            L" SuperNEMO @ 90% CL $$", 
            L" best $^{82}$Se",
            L" best world $$"
        ], 
        patchsize = (40, 40),
        L"SuperNEMO $0\nu\beta\beta \langle \lambda \rangle$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure_0nuRH_lambda.png"), f, px_per_unit = 5)
    f
    
end

### RH eta
### For eta: [25,43-38,81] 1E23 years

t_current_RH_e = get_sensitivities_vs_time(
    signal,
    background,
    SNparams;
    effFactor = 0.888
)

let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated sensitivity (yr)$$",yscale = log10, limits = (0, maximum(exp), 5e22, 9e24), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{$0\nu\beta\beta \langle \eta \rangle$ }_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 5e22, 9e24))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta \langle \eta \rangle$ }$ decay",)



    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3],0.7), linewidth = 1.5)

    p1 = lines!(ax, exp, t_current_RH_e, color = colors[2], linewidth = 5.5)
    lines!(ax2, t |> collect, t_current_RH_e, color = colors[2], linewidth = 5.5)
    
    h1 = hlines!(ax, [2.2e23], color = colors[2], linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)
    b1 = band!(ax, [0,35], [9.93e23], [38.81e23], color = (colors[2], 0.4), label = L"best world: $\langle \eta \rangle$")
        

    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, h1, b1], 
        [
            L" SuperNEMO @ 90% CL $$", 
            L" best $^{82}$Se",
            L" best world $$"
        ], 
        patchsize = (40, 40),
        L"SuperNEMO $0\nu\beta\beta \langle \eta \rangle$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure_0nuRH_eta.png"), f, px_per_unit = 5)
    f
    
end

### M1
signalM1 = get_process("bb0nuM1_foil_bulk", all_processes) |> first
set_nTotalSim!( signalM1, 0.1e8 )
set_bins!(signalM1, Bin_low:bin_width:Bin_high+bin_width)

t_current_nu0M1 = get_sensitivities_vs_time(
    signalM1,
    background,
    SNparams
)

let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated sensitivity (yr)$$",yscale = log10, limits = (0, maximum(exp), 1e22, 3e23), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{$0\nu\beta\beta \langle \eta \rangle$ }_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e22, 3e23))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta\chi^0$ }$ decay",)



    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3],0.7), linewidth = 1.5)

    p1 = lines!(ax, exp, t_current_nu0M1, color = colors[7], linewidth = 5.5)
    lines!(ax2, t |> collect, t_current_nu0M1, color = colors[7], linewidth = 5.5)
    
    h1 = hlines!(ax, [1.2e23], color = colors[7], linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)
    # b1 = band!(ax, [0,35], [9.93e23], [38.81e23], color = (colors[7], 0.4), label = L"best world: $\langle \eta \rangle$")
        

    text!(17.5/6.25, 2e22, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e22, 1e23, 1e24, 1e25], [L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e22, 1e23, 1e24, 1e25], [L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, h1], 
        [
            L" SuperNEMO @ 90% CL $$", 
            L" best $^{82}$Se",
        ], 
        patchsize = (40, 40),
        L"SuperNEMO $0\nu\beta\beta\chi^0$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure_0nuM1.png"), f, px_per_unit = 5)
    f
    
end

### M2

signalM2 = get_process("bb0nuM2_foil_bulk", all_processes) |> first
set_nTotalSim!( signalM2, 0.1e8 )
set_bins!(signalM2, Bin_low:bin_width:Bin_high+bin_width)

t_current_nu0M2 = get_sensitivities_vs_time(
    signalM2,
    background,
    SNparams
)

let exp = t .* SNparams["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated sensitivity (yr)$$",yscale = log10, limits = (0, maximum(exp), 1e21, 7e22), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{$0\nu\beta\beta \langle \eta \rangle$ }_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e21, 7e22))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta\chi^0\chi^0$ }$ decay",)



    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3], 0.7), linewidth = 1.5)

    p1 = lines!(ax, exp, t_current_nu0M2, color = colors[5], linewidth = 5.5)
    lines!(ax2, t |> collect, t_current_nu0M2, color = colors[5], linewidth = 5.5)
    
    h1 = hlines!(ax, [1.4e22], color = colors[5], linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)
    # b1 = band!(ax, [0,35], [9.93e23], [38.81e23], color = (colors[5], 0.4), label = L"best world: $\langle \eta \rangle$")
        

    text!(17.5/6.25, 2e22, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e21,1e22, 1e23, 1e24, 1e25], [L"10^{21}", L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e21,1e22, 1e23, 1e24, 1e25], [L"10^{21}", L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, h1], 
        [
            L" SuperNEMO @ 90% CL $$", 
            L" best $^{82}$Se",
        ], 
        patchsize = (40, 40),
        L"SuperNEMO $0\nu\beta\beta\chi^0$ sensitivity",
        )
    save(scriptsdir("LSM_report_2025_10_Christine", "nu0_sens_figure_0nuM2.png"), f, px_per_unit = 5)
    f
    
end

using DataFramesMeta, CSV
data_all = DataFrame(
    t = t |> collect,
    nu0bb = nu0_t_sens,
    nu0bb_RH_eta = t_current_RH_e,
    nu0bb_RH_lambda = t_current_RH_L,
    nu0bb_M1 = t_current_nu0M1,
    nu0bb_M2 = t_current_nu0M2
)

CSV.write("/home/maros/Work/Phd/SNSensitivityEstimate/scripts/LSM_report_2025_10_Christine/for_Yves/data.csv", data_all)