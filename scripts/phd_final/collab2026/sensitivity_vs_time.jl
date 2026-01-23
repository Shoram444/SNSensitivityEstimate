using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using ColorSchemes, SNSensitivityEstimate, CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings,StatsBase, FHist, Distributions
using FHist, ColorSchemes

colors = ["#041E42","#BE4D00","#951272","#006630","#005C8A","#FFB948","#605643","#302D23"]

df_nu0 = DataFrame(
    signal = [ L"0\nu\beta\beta", L"0\nu\beta\beta", L"0\nu\beta\beta", ],
    radon_level = [L"159 \muBq/m^3", L"2mBq/m^3", L"0.6 mBq/m^3" ],
    sens = [4.03e24, 2.88e24, 3.74e24],
    eff = [0.163, 0.163, 0.163],
    bkg = [1.46, 4.56, 1.9]
)

t = 0.1:0.1:5.0
α = 1.64485362695147
const sn_params = SNSensitivityEstimate.SNparams

function sensitivity_in_time(efficiency , background, time ;α = α, efffactor = 1.0, params = sn_params)
    get_tHalf(
        sn_params["W"],
        sn_params["foilMass"],
        sn_params["Nₐ"],
        time,
        sn_params["a"],
        efficiency * efffactor,
        (background)/ sn_params["tYear"] * time,
        α;
        approximate="table"
    )
end

sens_radon = [ map(x-> sensitivity_in_time(r.eff, r.bkg, x), t) for r in eachrow(df_nu0) ]


let exp = t .* sn_params["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e23, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{0\nu\beta\beta}$ decay",)


    h1 = hlines!(ax, [4.6e24], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    vlines!(ax, [6.25], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, sens_radon[1], color = colors[1], linewidth = 5.5)
    lines!(ax2, t |> collect, sens_radon[1], color = colors[1], linewidth = 5.5)
    
    # s1 = scatter!(ax, [17.5], [4.64e24], marker = :circle, markersize = 22, color = colors[2])
    
    p2 = lines!(ax, exp, sens_radon[2], color = colors[4], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, sens_radon[2], color = colors[4], linewidth = 5.5, linestyle = :dash)

    p3 = lines!(ax, exp, sens_radon[3], color = colors[5], linewidth = 5.5, linestyle = :dashdot)
    lines!(ax2, t |> collect, sens_radon[3], color = colors[5], linewidth = 5.5, linestyle = :dashdot)

    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    text!(1.0, 2e23, text = "6.25 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, p2, p3, h1], 
        [
            # L"$$  bayesian limit @ 90% CI",
            L"$$  Rn @ 0.15~$\textrm{mBq/m^3}$", 
            L"$$  Rn @ 2~$\textrm{mBq/m^3}$", 
            L"$$  Rn @ 0.6~$\textrm{mBq/m^3}$", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (60, 70),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity \\ frequentist limit @ 90% CL \\",
        )
    save(scriptsdir("phd_final/collab2026/figs","nu0_sens_figure.png"), f, px_per_unit = 5)
    f
    
end

df_nu0M1 = DataFrame(
    signal = [ L"0\nu\beta\beta\chi^2", L"0\nu\beta\beta\chi^2", L"0\nu\beta\beta\chi^2", ],
    radon_level = [L"159 \muBq/m^3", L"2mBq/m^3", L"0.6 mBq/m^3" ],
    sens = [1.56e23, 1.43e23, 1.52e23],
    eff = [0.019, 0.019, 0.019],
    bkg = [31.2, 37.7, 32.7]
)

sens_radonM1 = [ map(x-> sensitivity_in_time(r.eff, r.bkg, x), t) for r in eachrow(df_nu0M1) ]


let exp = t .* sn_params["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e22, 3e23), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, 5, 1e22, 3e23))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta\chi^0$ }$ decay",)


    h1 = hlines!(ax, [1.2e23], color = colors[7], linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    vlines!(ax, [6.25], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, sens_radonM1[1], color = colors[2], linewidth = 5.5)
    lines!(ax2, t |> collect, sens_radonM1[1], color = colors[2], linewidth = 5.5)
    
    # s1 = scatter!(ax, [17.5], [4.64e24], marker = :circle, markersize = 22, color = colors[2])
    
    p2 = lines!(ax, exp, sens_radonM1[2], color = colors[6], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, sens_radonM1[2], color = colors[6], linewidth = 5.5, linestyle = :dash)

    p3 = lines!(ax, exp, sens_radonM1[3], color = colors[7], linewidth = 5.5, linestyle = :dashdot)
    lines!(ax2, t |> collect, sens_radonM1[3], color = colors[7], linewidth = 5.5, linestyle = :dashdot)

    text!(17.5/6.25, 0.2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    text!(1.0, 0.2e23, text = "6.25 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e22, 1e23, 1e24, 1e25], [L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e22, 1e23, 1e24, 1e25], [L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, p2, p3, h1], 
        [
            # L"$$  bayesian limit @ 90% CI",
            L"$$  Rn @ 0.15~$\textrm{mBq/m^3}$", 
            L"$$  Rn @ 2~$\textrm{mBq/m^3}$", 
            L"$$  Rn @ 0.6~$\textrm{mBq/m^3}$", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (60, 70),
        L"SuperNEMO $$0\nu\beta\beta\chi^0$$ sensitivity \\ frequentist limit @ 90% CL",
        )
    save(scriptsdir("phd_final/collab2026/figs","nu0M1_sens_figure.png"), f, px_per_unit = 5)
    f
    
end


df_nu0M2 = DataFrame(
    signal = [ L"0\nu\beta\beta\chi^2\chi^2", L"0\nu\beta\beta\chi^2\chi^2", ],
    radon_level = [L"159 \muBq/m^3", L"2mBq/m^3" ],
    sens = [2.52e22, 2.50e22],
    eff = [0.101, 0.101],
    bkg = [47090, 47200]
)

sens_radonM2 = [ map(x-> sensitivity_in_time(r.eff, r.bkg, x), t) for r in eachrow(df_nu0M2) ]


let exp = t .* sn_params["foilMass"] |> collect
    f = Figure(size = (1000, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e21, 6e22), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, 5, 1e21, 6e22))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta\chi^0\chi^0$ }$ decay",)


    h1 = hlines!(ax, [1.4e22], color = (:black, 0.7), linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    vlines!(ax, [6.25], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, sens_radonM2[1], color = colors[6], linewidth = 5.5)
    lines!(ax2, t |> collect, sens_radonM2[1], color = colors[6], linewidth = 5.5)
    
    p2 = lines!(ax, exp, sens_radonM2[2], color = colors[4], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, sens_radonM2[2], color = colors[4], linewidth = 5.5, linestyle = :dash)


    text!(17.5/6.25, 0.6e22, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    text!(1.0, 0.6e22, text = "6.25 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e21, 1e22, 1e23, 1e24], [L"10^{21}",L"10^{22}", L"10^{23}", L"10^{24}"])
    ax2.yticks = ([1e21, 1e22, 1e23, 1e24], [L"10^{21}",L"10^{22}", L"10^{23}", L"10^{24}"])

    Legend(
        f[1,2], 
        [ p1, p2, h1], 
        [
            # L"$$  bayesian limit @ 90% CI",
            L"$$  Rn @ 0.15~$\textrm{mBq/m^3}$", 
            L"$$  Rn @ 2~$\textrm{mBq/m^3}$", 
            # L"$$  Rn @ 0.6~$\textrm{mBq/m^3}$", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (60, 70),
        L"SuperNEMO $$0\nu\beta\beta\chi^0\chi^0$$ sensitivity \\ frequentist limit @ 90% CL",
        )
    save(scriptsdir("phd_final/collab2026/figs","nu0M2_sens_figure.png"), f, px_per_unit = 5)
    f
    
end




df_nu0_nd = DataFrame(
    signal = [ L"0\nu\beta\beta", L"0\nu\beta\beta"],
    radon_level = [L"150 \muBq/m^3", L"2mBq/m^3"],
    sens = [4.24e24, 3.51e24],
    eff = [0.19, 0.19],
    bkg = [1.81, 3.71]
)

sens_radon_nd = [ map(x-> sensitivity_in_time(r.eff, r.bkg, x), t) for r in eachrow(df_nu0_nd) ]


let exp = t .* sn_params["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 1e23, 1e25))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{0\nu\beta\beta}$ decay",)


    h1 = hlines!(ax, [4.6e24], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)

    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    vlines!(ax, [6.25], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, sens_radon[1], color = colors[1], linewidth = 5.5)
    lines!(ax2, t |> collect, sens_radon[1], color = colors[1], linewidth = 5.5)
    
    # s1 = scatter!(ax, [17.5], [4.64e24], marker = :circle, markersize = 22, color = colors[2])
    
    p2 = lines!(ax, exp, sens_radon[2], color = colors[4], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, sens_radon[2], color = colors[4], linewidth = 5.5, linestyle = :dash)

    p3 = lines!(ax, exp, sens_radon_nd[1], color = colors[2], linewidth = 5.5, )
    lines!(ax2, t |> collect, sens_radon_nd[1], color = colors[2], linewidth = 5.5, )

    p4 = lines!(ax, exp, sens_radon_nd[2], color = colors[3], linewidth = 5.5, linestyle = :dash)
    lines!(ax2, t |> collect, sens_radon_nd[2], color = colors[3], linewidth = 5.5, linestyle = :dash)


    text!(17.5/6.25, 2e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    text!(1.0, 2e23, text = "6.25 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1, p2, p3, p4, h1], 
        [
            # L"$$  bayesian limit @ 90% CI",
            L"$$ 1D ROI Rn @ 0.15~$\textrm{mBq/m^3}$", 
            L"$$ 1D ROI Rn @ 2~$\textrm{mBq/m^3}$", 
            L"$$ ND ROI Rn @ 0.15~$\textrm{mBq/m^3}$", 
            L"$$ ND ROI Rn @ 2~$\textrm{mBq/m^3}$", 
            L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (60, 70),
        L"SuperNEMO $$0\nu\beta\beta$$ sensitivity \\ frequentist limit @ 90% CL \\",
        )
    save(scriptsdir("phd_final/collab2026/figs","nd_nu0_sens_figure.png"), f, px_per_unit = 5)
    f
    
end



### RH
df_RH = DataFrame(
    signal = [ L"\nu_L\nu_R\beta\beta" ],
    radon_level = [L"150 \muBq/m^3",],
    sens = [1.41e22, ],
    eff = [0.103],
    bkg = [156000]
)

sens_RH = [ map(x-> sensitivity_in_time(r.eff, r.bkg, x), t) for r in eachrow(df_RH) ]



let exp = t .* sn_params["foilMass"] |> collect
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity (yr)",yscale = log10, limits = (0, maximum(exp), 1e21, 0.5e23), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    # ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated $T^{0\nu\beta\beta}_{1/2}$ sensitivity @ 90% CL (yr)",yscale = log10, limits = (0, maximum(exp), 1e23, 1e25), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, 5, 1e21, 0.5e23))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$\nu_L\nu_R\beta\beta$ }$ decay",)



    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    vlines!(ax, [6.25], linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
    
    p1 = lines!(ax, exp, sens_RH[1], color = :black, linewidth = 5.5)
    lines!(ax2, t |> collect, sens_RH[1], color = :black, linewidth = 5.5)
    

    text!(17.5/6.25, 0.2e22, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    text!(1.0, 0.2e22, text = "6.25 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))
    
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = ([1e21, 1e22, 1e23, 1e24, 1e25], [L"10^{21}",L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e21, 1e22, 1e23, 1e24, 1e25], [L"10^{21}",L"10^{22}",L"10^{23}", L"10^{24}", L"10^{25}"])

    Legend(
        f[1,2], 
        [ p1], 
        [
            # L"$$  bayesian limit @ 90% CI",
            L"$$  Rn @ 0.15~$\textrm{mBq/m^3}$", 
            # L"$$  Rn @ 2~$\textrm{mBq/m^3}$", 
            # L"$$  Rn @ 0.6~$\textrm{mBq/m^3}$", 
            # L"$$  CUPID-0 @ 90% CI\n bayesian limit"
        ], 
        patchsize = (60, 70),
        L"SuperNEMO $\nu_L\nu_R\beta\beta$ sensitivity \\ frequentist limit @ 90% CL",
        )
    save(scriptsdir("phd_final/collab2026/figs","RH_sens_figure.png"), f, px_per_unit = 5)
    f
    
end