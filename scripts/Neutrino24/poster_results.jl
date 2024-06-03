using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Distributions
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

processes = load_processes("fal5_12perc", "sumE")

# scaled to Bq

names = [
    "bb_foil_bulk",
    "Bi214_foil_bulk",
    "Tl208_foil_bulk",
    # "K40_foil_bulk",
    # "Pa234m_foil_bulk",
    "Bi214_wire_surface",
    # "Bi214_hall_surface",
    # "Tl208_hall_surface",
    # "K40_hall_surface",
]

bkgActivitiesDict= Dict(
    :bb_foil_bulk => SigActivityParams[:bb_foil_bulk].val * SNparams["foilMass"],
    :Bi214_foil_bulk => BkgActivityParams[:Bi214_foil_bulk] * SNparams["foilMass"],
    :Tl208_foil_bulk => BkgActivityParams[:Tl208_foil_bulk] * SNparams["foilMass"],
    :K40_foil_bulk => BkgActivityParams[:K40_foil_bulk].val * SNparams["foilMass"],
    :Pa234m_foil_bulk => BkgActivityParams[:Pa234m_foil_bulk].val * SNparams["foilMass"],
    :Bi214_wire_surface => BkgActivityParams[:Bi214_wire_surface] * SNparams["gasVolume"],
    :Bi214_hall_surface => BkgActivityParams[:Bi214_hall_surface] * SNparams["hall_surface"],
    :Tl208_hall_surface => BkgActivityParams[:Tl208_hall_surface] * SNparams["hall_surface"],
    :K40_hall_surface => BkgActivityParams[:K40_hall_surface] * SNparams["hall_surface"],
)

nTotalSimDict = Dict(
    :bb_foil_bulk => 267e6,
    :Bi214_foil_bulk => 3*99e6,
    :Tl208_foil_bulk => 3*98e6,
    :K40_foil_bulk => 0,
    :Pa234m_foil_bulk => 0,
    :Bi214_wire_surface => 3*79e6,
    :Bi214_hall_surface => 0,
    :Tl208_hall_surface => 0,
    :K40_hall_surface => 0,
)

bw = 50

h1d_bb_foil_bulk = get_estimated_counts("bb_foil_bulk", bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
Bi214_foil_bulk = get_estimated_counts("Bi214_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
Tl208_foil_bulk = get_estimated_counts("Tl208_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
# K40_foil_bulk = get_estimated_counts("K40_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
# Pa234m_foil_bulk = get_estimated_counts("Pa234m_foil_bulk" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
Bi214_wire_surface = get_estimated_counts("Bi214_wire_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
# Bi214_hall_surface = get_estimated_counts("Bi214_hall_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
# Tl208_hall_surface = get_estimated_counts("Tl208_hall_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)
# K40_hall_surface = get_estimated_counts("K40_hall_surface" , bkgActivitiesDict, nTotalSimDict, SNparams["t"];  bins = 0:bw:3500)

bkg_hists = [
    h1d_bb_foil_bulk,
    Bi214_foil_bulk + Tl208_foil_bulk,  # + K40_foil_bulk + Pa234m_foil_bulk
    Bi214_wire_surface,
    Tl208_foil_bulk
    # Bi214_hall_surface + Tl208_hall_surface + K40_hall_surface
]

min_cts = minimum(@. minimum( filter(x-> x>0, bincounts(bkg_hists)) ) );
colors = ["#003865", "#FFB948", "#52473B", "#9A3A06", ]

let
	f = Figure(size = (800, 600), fontsize= 26, fonts = (; regular = "TeX"))
	ax = Axis(f[1:2,1], xlabel = "electron sum energy [keV]", ylabel = "counts / $bw keV" , title = "simulated SuperNEMO background", yscale =log10, xticklabelrotation=45)
	ax2 = Axis(f[2,2], xlabel = "energy [keV]", title = L"\textrm{0\nu\beta\beta ROI}", yscale =log10, limits=(2700,3000, nothing,nothing), yaxisposition = :right, xticklabelrotation=45)
	
	labels= [L"2\nu\beta\beta", "internal", "radon", "external"]
	st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
	for i=2:length(bkg_hists)
		hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
	end

	st2 = hist!(ax2, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
	for i=2:length(bkg_hists)
		hist!(ax2, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
	end
	
    text!(ax, 2000, 1.5e1, text= "Preliminary", fontsize = 40, rotation = pi/6, align = (:center, :baseline), color = (:red, 0.6))
	ylims!(ax, min_cts, 1e6)
	xlims!(ax, 200, 3200)
	ylims!(ax2, min_cts, 1e1)
    ax.xticks = (collect(500:500:3000), string.(collect(500:500:3000)))
    ax2.xticks = ([2700, 3000], ["2700", "3000"])
	f[1, 2] = Legend(f, ax)
	colgap!(f.layout, 1, Relative(0.05))
    save(scriptsdir("Neutrino24","fig.png"), f)
	f
end