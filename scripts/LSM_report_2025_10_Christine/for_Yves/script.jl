using CSV, DataFramesMeta, LaTeXStrings

data = CSV.File("/home/maros/Work/Phd/SNSensitivityEstimate/scripts/LSM_report_2025_10_Christine/for_Yves/data.csv") |> DataFrame
t = data.t
exposure = t .* 6.25

nu0bb_RH_eta = data.nu0bb_RH_eta
nu0bb_RH_lambda = data.nu0bb_RH_lambda

colors = ["#003865", "#FFB948", "#52473B", "#9A3A06", "#951272", "#006630", "#005C8A", "#FFB948"]

# figure for eta
let 
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated sensitivity (yr)$$",yscale = log10, limits = (0, maximum(exposure), 5e22, 9e24), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 5e22, 9e24))

    # top label
    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta \langle \eta \rangle$ }$ decay",)

    # line for the 17.5 kg.yr exposure
    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3],0.7), linewidth = 1.5)
    text!(17.5/6.25, 1e23, text = "17.5 kg.yr", align = (:center, :top), rotation = pi/2, fontsize = 20, color = (:black, 0.3))

    # line for the sensitivity
    p1 = lines!(ax, exposure, nu0bb_RH_eta, color = colors[2], linewidth = 5.5)
    lines!(ax2, t |> collect, nu0bb_RH_eta, color = colors[2], linewidth = 5.5)
    
    # best limits
    h1 = hlines!(ax, [2.2e23], color = colors[2], linestyle = :dash, label = L"best $^{82}$Se: $\langle \eta \rangle$", linewidth = 2)
    b1 = band!(ax, [0,35], [9.93e23], [38.81e23], color = (colors[2], 0.4), label = L"best world: $\langle \eta \rangle$")
        
    # hide unneeded axis decorations
    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    # manually set ticks
    ax.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])
    ax2.yticks = ([1e23, 1e24, 1e25], [L"10^{23}", L"10^{24}", L"10^{25}"])

    # manually set legend
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

    f    
end


# figure for lambda
let 
    f = Figure(size = (950, 600), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
    ax = Axis(f[1, 1], xlabel = L"$$ exposure (kg.yr)",  ylabel = L"Simulated sensitivity (yr)$$",yscale = log10, limits = (0, maximum(exposure), 5e22, 5e24), yminorticks=IntervalsBetween(10), yminorticksvisible = true, yminorticksize = 4)
    ax2 = Axis(f[1, 1], xlabel = L"$$time (yr)", yscale = log10, xaxisposition = :top, limits = (0, maximum(t), 5e22, 5e24))

    Label(f[0,1], text = L"Simulated SuperNEMO sensitivity to $\textrm{$0\nu\beta\beta \langle \lambda \rangle$ }$ decay",)


    vlines!(ax, [17.5], linestyle = (:dash, :dense), color = (colors[3],0.7), linewidth = 1.5)

    p1 = lines!(ax, exposure, nu0bb_RH_lambda, color = colors[6], linewidth = 5.5)
    lines!(ax2, t |> collect, nu0bb_RH_lambda, color = colors[6], linewidth = 5.5)
    
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
    f
    
end