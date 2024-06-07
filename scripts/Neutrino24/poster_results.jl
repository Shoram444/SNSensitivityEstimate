using DrWatson

@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie, UnROOT, ColorSchemes,DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Distributions
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

function get_estimated_counts(data::Vector{<:Real}, activity::Real, nTotalSim::Real, time::Real; bins = 0:50:3500)
    h = Hist1D(; binedges = bins, counttype=Float64)
    efficiency = length(data) / nTotalSim
    w = activity*efficiency*time 
    for d in data
        push!(h, d)
    end
    h = normalize(h, width =false)

    h *= w
    return h
end

function get_estimated_counts(isotopeName::String, bkgActivitiesDict::Dict, nTotalSimDict::Dict, time::Real;  bins = 0:50:3500)
    data = get_process(isotopeName, processes).dataVector
    return get_estimated_counts(data, bkgActivitiesDict[Symbol(isotopeName)], nTotalSimDict[Symbol(isotopeName)], time; bins = 0:50:3500)
end

processes = load_processes("fal5_12perc_Boff", "sumE")

# scaled to Bq

names = [
    "bb_foil_bulk",
    "Bi214_foil_bulk",
    "Tl208_foil_bulk",
    "K40_foil_bulk",
    "Pa234m_foil_bulk",
    "Bi214_wire_surface",
    "Bi214_hall_surface",
    "Tl208_hall_surface",
    "K40_hall_surface",
]

bkgActivitiesDict= Dict(
    :bb_foil_bulk => SigActivityParams[:bb_foil_bulk].val * SNparams["foilMass"],
    :Bi214_foil_bulk => BkgActivityParams[:Bi214_foil_bulk] * SNparams["foilMass"],
    :Tl208_foil_bulk => BkgActivityParams[:Tl208_foil_bulk] * SNparams["foilMass"],
    :K40_foil_bulk => BkgActivityParams[:K40_foil_bulk] * SNparams["foilMass"],
    :Pa234m_foil_bulk => BkgActivityParams[:Pa234m_foil_bulk] * SNparams["foilMass"],
    :Bi214_wire_surface => BkgActivityParams[:Bi214_wire_surface] * SNparams["gasVolume"],
    :Bi214_hall_surface => BkgActivityParams[:Bi214_hall_surface] * SNparams["hall_surface"],
    :Tl208_hall_surface => BkgActivityParams[:Tl208_hall_surface] * SNparams["hall_surface"],
    :K40_hall_surface => BkgActivityParams[:K40_hall_surface] * SNparams["hall_surface"],
)

nTotalSimDict = Dict(
    :bb_foil_bulk => 267e6,
    :Bi214_foil_bulk => 3*99e6,
    :Tl208_foil_bulk => 3*98e6,
    :K40_foil_bulk => 958*3e5,
    :Pa234m_foil_bulk => 3e8,
    :Bi214_wire_surface => 3*79e6,
    :Bi214_hall_surface => 1e9,
    :Tl208_hall_surface => 1e9,
    :K40_hall_surface => 1e9,
)

bw = 50

h1d_bb_foil_bulk = get_estimated_counts("bb_foil_bulk", bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_Bi214_foil_bulk = get_estimated_counts("Bi214_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_Tl208_foil_bulk = get_estimated_counts("Tl208_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_K40_foil_bulk = get_estimated_counts("K40_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_Pa234m_foil_bulk = get_estimated_counts("Pa234m_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_Bi214_wire_surface = get_estimated_counts("Bi214_wire_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_Bi214_hall_surface = get_estimated_counts("Bi214_hall_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_Tl208_hall_surface = get_estimated_counts("Tl208_hall_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
h1d_K40_hall_surface = get_estimated_counts("K40_hall_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)

bkg_hists = [
    h1d_bb_foil_bulk,
    h1d_Bi214_foil_bulk + h1d_Tl208_foil_bulk + h1d_K40_foil_bulk + h1d_Pa234m_foil_bulk,
    h1d_Bi214_wire_surface,
    h1d_Bi214_hall_surface + h1d_Tl208_hall_surface + h1d_K40_hall_surface + h1d_Tl208_foil_bulk*inv(1000)
]

min_cts = minimum(@. minimum( filter(x-> x>0, bincounts(bkg_hists)) ) ) *10
colors = colorschemes[:julia] #["#003865", "#FFB948", "#52473B", "#9A3A06", ]
# colors = ["#003865", "#FFB948", "#52473B", "#9A3A06", ]

with_theme(theme_latexfonts()) do 
	f = Figure(size = (800, 600), fontsize = 24, fonts = (; regular = "TeX"), figure_padding = 24)
	ax = Axis(f[1:2,1], xlabel = "Summed 2-electron energy [keV]", ylabel = "Counts per $bw keV" , title = "Simulated SuperNEMO background; 17.5 kg.yr", yscale =log10, xticklabelrotation=45)
	ax2 = Axis(f, bbox=BBox(590,745,410,530), yscale =log10,yaxisposition = :left,xticklabelrotation=0, xticklabelsize = 20, backgroundcolor=(:blue, 0.1)) #, title = L"\textrm{0\nu\beta\beta ROI}"
	
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
    lines!(ax, [2700, 2700, 2555], [min_cts, 1e1, 0.45e2], color=(:black), linewidth=3.5)
    lines!(ax, [3200, 3200, 3340], [min_cts, 1e1, 0.45e2], color=(:black), linewidth=3.5)
    text!(ax, 2600, 7e4, text=L"\textrm{0\nu\beta\beta ROI}", fontsize=32)

    scatter!(ax, [SNparams["Q"]], [min_cts*1.9], marker= :dtriangle, markersize = 18, color=:red)

    text!(ax, 1900, 1.5e-1, text= "Preliminary", fontsize = 40, rotation = pi/6, align = (:center, :baseline), color = (:red, 0.8))
	ylims!(ax, min_cts, 1e6)
	xlims!(ax, 200, 3500)
	ylims!(ax2, 1e-4, 5e2)
	xlims!(ax2, 2700, 3200)
    ax.xticks = (collect(500:500:3900), string.(collect(500:500:3900)))
    ax.yticks = ([1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
    ax2.xticks= ( [2700, 3200], ["2700", "3200"] )
    ax2.yticks= ( [1e-4, 1e-2, 1e0, 1e2], [L"10^{-4}", L"10^{-2}", L"10^{0}", L"10^{2}"] )
    
    axislegend(ax, position = :lt, orientation = :horizontal, nbanks=2)
    
    save(scriptsdir("Neutrino24","fig.png"), f)
    save(scriptsdir("Neutrino24","fig.svg"), f)
	f
end




signal = get_process("bb0nu_foil_bulk", processes)
background = [
    get_process("bb_foil_bulk", processes),
    get_process("Bi214_foil_bulk", processes),
    get_process("Bi214_wire_surface", processes),
    get_process("Tl208_foil_bulk", processes),
    get_process("Bi214_hall_surface", processes),
    get_process("Tl208_hall_surface", processes),
    get_process("K40_hall_surface", processes),
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)
set_nTotalSim!( signal, 99e6 )
set_nTotalSim!( background[1], 267e6 )
set_nTotalSim!( background[2], 3*99e6 )
set_nTotalSim!( background[3], 3*79e6 )
set_nTotalSim!( background[4], 3*98e6 )

Q_keV = SNparams["Q"]
α = 1.64485362695147


t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)