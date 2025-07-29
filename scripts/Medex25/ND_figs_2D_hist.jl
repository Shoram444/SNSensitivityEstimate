using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, FHist, ColorSchemes
colors          = ["#003865", "#FFB948", "#52473B", "#9A3A06", ]


include(scriptsdir("ND/results/best_rois.jl"))
include(srcdir("params/Params.jl"))

vars = [
    "phi", 
    "sumE",
    "singleE",
    "lPint",
    "lPext",
    "dy",
    "dz",
]
bins = (
    phi = (0,180),
    sumE = (0, 3500),
    singleE = (0, 3500),
    lPint = (0,50),
    lPext = (0,50),
    dy = (0,250),
    dz = (0,250),
)



processesND = load_ndim_processes(datadir("mva/fal5_TKrec_J40"), bins, vars)


signalNDM1 = get_process("bb0nu_foil_bulk", processesND) |> first
set_nTotalSim!( signalNDM1, 0.1e8 )
phiSignalNDM1 = getproperty.(signalNDM1.data,:phi)
esumSignalNDM1 = getproperty.(signalNDM1.data,:sumE)

backgroundNDbb = get_process("bb_foil_bulk", processesND) |> first
set_signal!(backgroundNDbb, false)
backgroundND = get_process("Tl208_foil_bulk", processesND) |> first

phiBackgroundND = vcat(getproperty.(backgroundND.data,:phi), getproperty.(backgroundNDbb.data,:phi))
esumBackgroundND = vcat(getproperty.(backgroundND.data,:sumE), getproperty.(backgroundNDbb.data,:sumE))

h2dsignal = Hist2D(
        (phiSignalNDM1, esumSignalNDM1);
        binedges = (0:5:180, 0:100:3500)
    ) 
h2dsignal.bincounts .*= inv(0.1e8)
h2dbackground = Hist2D(
        (phiBackgroundND, esumBackgroundND);
        binedges = (0:5:180, 0:100:3500)
    ) |> normalize
h2dbackground.bincounts .*= 1e5

min_cts_sig = minimum(filter(x -> x > 0, h2dsignal.bincounts))
min_cts_bkg = minimum(filter(x -> x > 0, h2dbackground.bincounts))

with_theme(theme_latexfonts()) do 
	f = Figure(size = (1000, 400), fontsize = 22, fonts = (; regular = "TeX"), figure_padding = 20)
	
	a1 = Axis(f[1,1], xlabel = L"$\varphi$", ylabel = L"$E_{sum}$", title = L"Signal efficiency $$")
	

    c = Makie.to_colormap([colors[1],colors[2]])
    cmap_sig = cgrad(c, 0.1*min_cts_sig, scale = :log10)
    cmap_bkg = cgrad(c, 0.1*min_cts_bkg, scale = :log10)

	h1 = plot!(a1, h2dsignal, colormap = cmap_sig)
	Colorbar(f[1,2], h1,)

	a2 = Axis(f[1,3], xlabel = L"$\varphi$", ylabel = L"$E_{sum}$", title = L"Background counts $$")
	h2 = plot!(a2, h2dbackground, colormap = cmap_bkg)
	Colorbar(f[1,4], h2,)


	poly!(
		a1,
		[10,175,175,10],
		[2500,2500,3100,3100],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
    )
    text!(a1, 17,2800, text = "2D - ROI", align = (:left, :top), fontsize = 20, color = (:red, 1))
    
	poly!(
		a2,
		[10,175,175,10],
		[2500,2500,3100,3100],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
    )
    text!(a2, 17,2800, text = "2D - ROI", align = (:left, :top), fontsize = 20, color = (:red, 1))
	
	save(scriptsdir("Medex25", "ND_phi_sumE.png"), f, px_per_unit = 5)
	f
end