using DrWatson
@quickactivate "SNSensitivityEstimate"

using SNSensitivityEstimate, CairoMakie, FHist, ColorSchemes


include(scriptsdir("ND/results/best_rois.jl"))
include(srcdir("params/Params.jl"))

analysisDict    = Dict(:signal => "bb0nu_foil_bulk", :bining => (300, 3500), :bin_width => 100, :mode => "sumE")

roi             = bb0nu_roi
roi[:sumE]      = (300, 3500)

Bin_low, Bin_high, bin_width = analysisDict[:bining][1], analysisDict[:bining][2], analysisDict[:bin_width]


all_processes   = load_data_processes(datadir("mva/fal5_TKrec_J40"), analysisDict[:mode], fwhm=0.0, roi = topology_roi)
signal          = get_process(analysisDict[:signal], all_processes) |> first
background      = [get_process(b, all_processes) |> first for b in backgrounds]


for b in background
    set_bins!(b, Bin_low:bin_width:Bin_high+bin_width)
end
set_bins!(signal, Bin_low:bin_width:Bin_high+bin_width)


set_signal!(background[1], false)

bkg_hist        = get_bkg_counts_1D.(background)
bkg_hists       = [bkg_hist[1], sum(vcat(bkg_hist[2], bkg_hist[4:6])), bkg_hist[3], bkg_hist[7]]


min_cts         = minimum(@. minimum( filter(x-> x>0, bincounts(bkg_hists)) ) ) 
colors          = ["#003865", "#FFB948", "#52473B", "#9A3A06", ]

bw              = analysisDict[:bin_width]

# logscale
with_theme(theme_latexfonts()) do 
	f = Figure(size = (800, 600), fontsize = 24, fonts = (; regular = "TeX"), figure_padding = 24)
	ax = Axis(f[1:2,1], xlabel = "Summed 2-electron energy (keV)", ylabel = "Counts per $bw keV" , title = "Simulated SuperNEMO background; 17.5 kg.yr", yscale =log10, xticklabelrotation=45)
	ax2 = Axis(f, bbox=BBox(550,745,410,530), yscale =log10,yaxisposition = :left,xticklabelrotation=0, xticklabelsize = 20, backgroundcolor=(:blue, 0.1)) #, title = L"\textrm{0\nu\beta\beta ROI}"
	
    labels= [L"2\nu\beta\beta", "internal", "radon", "external"]
	st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
	for i=2:length(bkg_hists)
		hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
	end

	st2 = hist!(ax2, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
	for i=2:length(bkg_hists)
		hist!(ax2, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
	end
	
    lines!(ax, [SNparams["Q"], SNparams["Q"]], [min_cts*1.9, 0.12], color=(:red, 1), linewidth=3.5)
    text!(ax, SNparams["Q"], 0.13, text= L"\mathbf{\textrm{Q_{\beta\beta}}}", fontsize = 26, align = (:center, :baseline), color = (:red, 1))
    lines!(ax, [2700, 2700, 2455], [min_cts, 1e2, 0.6e3], color=(:black), linewidth=3.5)
    lines!(ax, [3100, 3100, 3340], [min_cts, 1e2, 0.6e3], color=(:black), linewidth=3.5)
    text!(ax, 2600, 2e5, text=L"\textrm{0\nu\beta\beta ROI}", fontsize=32)

    scatter!(ax, [SNparams["Q"]], [min_cts*1.9], marker= :dtriangle, markersize = 18, color=:red)

    # text!(ax, 1900, 1.5e-1, text= "Preliminary", fontsize = 40, rotation = pi/6, align = (:center, :baseline), color = (:red, 0.8))
	ylims!(ax, min_cts, 5e6)
	xlims!(ax, 300, 3500)
	ylims!(ax2, min_cts, 5e2)
	xlims!(ax2, 2700, 3100)
    ax.xticks = (collect(500:500:3900), string.(collect(500:500:3900)))
    ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
    ax2.xticks= ( [2700, 3100], ["2700", "3100"] )
    ax2.yticks= ( [1e-5, 1e-3, 1e-1, 1e1], [L"10^{-5}", L"10^{-3}", L"10^{-1}", L"10^{1}"] )
    
    axislegend(ax, position = :lt, orientation = :horizontal, nbanks=2)
    
    save(scriptsdir("Medex25","background_log.png"), f)
    save(scriptsdir("Medex25","background_log.svg"), f)
	f
end


max_cts = 1.7* maximum(@. maximum( filter(x-> x>0, bincounts(bkg_hists)) ) )
# linear scale
with_theme(theme_latexfonts()) do 
    f = Figure(size = (800, 600), fontsize = 24, fonts = (; regular = "TeX"), figure_padding = 30)
    ax = Axis(f[1:2,1], xlabel = "Summed 2-electron energy (keV)", ylabel = "Counts per $bw keV" , title = "Simulated SuperNEMO background; 17.5 kg.yr", xticklabelrotation=0, backgroundcolor = (:white, 1.0))
    ax2 = Axis(f, bbox=BBox(500,745,290,460), yaxisposition = :left,xticklabelrotation=0, xticklabelsize = 20, backgroundcolor = (:white, 1.0))
    
    labels= [L"2\nu\beta\beta", "internal", "radon", "external"]
    st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
    for i=2:length(bkg_hists)
        hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
    end

    st2 = hist!(ax2, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
    for i=2:length(bkg_hists)
        hist!(ax2, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
    end
    
    lines!(ax, [SNparams["Q"], SNparams["Q"]], [min_cts*1.9, 1600], color=(:red, 1), linewidth=3.5)
    text!(ax, SNparams["Q"], 1650, text= L"\mathbf{\textrm{Q_{\beta\beta}}}", fontsize = 26, align = (:center, :baseline), color = (:red, 1))
    lines!(ax, [2700,  2205], [min_cts*0.9,  6800], color=(:black), linewidth=3.5)
    lines!(ax, [3100,  3340], [min_cts*0.9,  6800], color=(:black), linewidth=3.5)
    text!(ax2, 2810,1.6, text=L"\textrm{0\nu\beta\beta ROI}", fontsize=28)
    # text!(ax, 2150, 6800, text="log-scale", fontsize=24)
    scatter!(ax, [SNparams["Q"]], [1.5e2], marker= :dtriangle, markersize = 18, color=:red)
    # text!(ax, 1900, 1.5e-1, text= "Preliminary", fontsize = 40, rotation = pi/6, align = (:center, :baseline), color = (:red, 0.8))
    ylims!(ax, min_cts, max_cts)
    xlims!(ax, 300, 3500)
    ylims!(ax2, 0, 2.3)#0.3e-1, 5e2)
    xlims!(ax2, 2700, 3100)
    ax.xticks = (collect(500:500:3900), string.(collect(500:500:3900)))
    # ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])       
    ax2.xticks= ( [2700, 3100], ["2700", "3100"] )
    ax2.yticks= ( [0, 1, 2], [L"0", L"1", L"2"] )
    axislegend(ax, position = :lt, orientation = :horizontal, nbanks=2)

    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax2, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    save(scriptsdir("Medex25","background_linear.png"), f)
    save(scriptsdir("Medex25","background_linear.svg"), f)
    f
end

