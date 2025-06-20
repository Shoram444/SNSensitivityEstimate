using DrWatson
@quickactivate "SNSensitivityEstimate"

using SNSensitivityEstimate, CairoMakie, FHist, ColorSchemes

colors          = ["#003865", "#FFB948", "#52473B", "#9A3A06", ]

nu0_ND = 4.11e24
nu0_1D = 3.53e24

nu0M1_ND = 1.51e23
nu0M1_1D = 1.30e23

nu0M2_1D = 2.55e22
nu0M2_ND = 2.9e22


with_theme(theme_latexfonts()) do 
    f = Figure(size = (1050, 500), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(
        f[1, 1], 
        ylabel = L"improvement (%) $$", 
        title = "Comparison of the sensitivities obtained from \nthe traditional 1D analysis vs the ND analysis",
        limits = (0, 4, -1, 21), 
        yminorticks=IntervalsBetween(10), 
        yminorticksvisible = true, 
        yminorticksize = 4, 
        xticklabelrotation = 0
        )

    snu01 = scatter!(ax, [1], [0],  color = colors[1], markersize = 35, marker = :circle)
    snu0N = scatter!(ax, [1], [nu0_ND/nu0_1D * 100 - 100], color = colors[1], markersize = 35, marker = :star5)
    snu0M1 = scatter!(ax, [2], [0], color = colors[2], markersize = 35, marker = :circle)
    snu0M1 = scatter!(ax, [2], [nu0M1_ND/nu0M1_1D * 100 - 100], color = colors[2], markersize = 35, marker = :star5)
    snu0M2 = scatter!(ax, [3], [0],  color = colors[3], markersize = 35, marker = :circle)
    snu0M2 = scatter!(ax, [3], [nu0M2_ND/nu0M2_1D * 100 - 100],  color = colors[3], markersize = 35, marker = :star5)

    p1 = scatter!(ax, [5], [0], color = (:black, 0.5), markersize = 35, marker = :circle)
    p2 = scatter!(ax, [5], [0], color = (:black, 0.5), markersize = 35, marker = :star5)
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.xticks = (1:3, [L"$0\nu\beta\beta$   ", L"$0\nu\beta\beta\chi^0$, n=1   ", L"$0\nu\beta\beta\chi^0\chi^0$, n=3   "])


    Legend(f[1,2], [p1,p2], ["1D analysis", "ND analysis"], patchsize = (40, 40))
    save(scriptsdir("Medex25", "sens_1D_vs_ND.png"), f, px_per_unit = 5)
    f
end

