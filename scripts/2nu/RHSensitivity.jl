using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Distributions
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))


## Phi
all_processes = load_processes("fal5_8perc", "phi")

signal_rfile = ROOTFile(datadir("sims","SDBDRC_vertex_prob","RH037_foil_bulk_SDBDRC_vertex_prob.root"))
signal_data = LazyTree(signal_rfile, "tree", ["phi"]) |> DataFrame


signal = Process(
    signal_data.phi,
    "RH_foil_bulk",
    true,    
    SigActivityParams[:RH037_foil_bulk],
    9.019296e7,
    1e7,
    0:5:180,
    6.07,
    )
background = [
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
]

set_signal!(background[1], false)

α = 1.64485362695147
t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)


h_b = get_bkg_counts_1D(background...)
h_bs = get_bkg_counts_1D.(background)
h_s = get_bkg_counts_1D(signal) * inv(100)

let 
    fig = Figure()
    a = Axis(fig[1,1], )
    labels = vcat(  "RH", [p.isotopeName for p in background], )

    all = vcat(h_s, h_bs, )
    st = hist!(a, sum(all), label =labels[1],c=1)
	for i=2:length(all)
		hist!(a, sum(all[i:end]), label=labels[i], c=i)
	end
    elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
    title = "Processes"
    Legend(fig[1,2], elements, labels, title)

    fig
end


let 
    fig = Figure()
    a = Axis(fig[1,1], yscale =log10, limits=(nothing, nothing, 1e-5, 1e4))
    labels = vcat(  "RH", [p.isotopeName for p in background], )

    all = vcat(h_s, h_bs, )
    st = hist!(a, sum(all), label =labels[1],c=1)
	for i=2:length(all)
		hist!(a, sum(all[i:end]), label=labels[i], c=i)
	end
    elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
    title = "Processes"
    Legend(fig[1,2], elements, labels, title)

    fig
end

let 
    fig = Figure()
    a = Axis(fig[1,1], )
    labels = vcat(  "RH", [p.isotopeName for p in background], )

    d_b = [p.dataVector for p in background]
    d_s = signal.dataVector
    p = hist!(a, d_s, color = (Makie.wong_colors()[1], 0.4), normalization = :pdf, linewidth = 4, scale_to=1.0)
	for i=2:length(d_b)+1
		stephist!(a, d_b[i-1], label=labels[i], color=Makie.wong_colors()[i], normalization = :none, linewidth = 4, scale_to=1.0)
	end
    elements = [PolyElement(polycolor = Makie.wong_colors()[i]) for i in 1:length(labels)]
    title = "Processes"
    Legend(fig[1,2], elements, labels, title)

    fig
end

plot(t12MapESum)

## Single
all_processes = load_processes("fal5_8perc", "singleE")

signal_rfile = ROOTFile(datadir("sims","SDBDRC_vertex_prob","RH037_foil_bulk_SDBDRC_vertex_prob.root"))
signal_data = LazyTree(signal_rfile, "tree", ["reconstructedEnergy1", "reconstructedEnergy2"]) |> DataFrame


signal = Process(
    vcat(signal_data.reconstructedEnergy1, signal_data.reconstructedEnergy2),
    "RH_foil_bulk",
    true,    
    1.0,
    9.019296e7,
    1e8,
    0:50:3500,
    6.07,
    )
background = [
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
]

set_signal!(background[1], false)

α = 1.64485362695147
t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)


## Sum
all_processes = load_processes("fal5_8perc", "sumE")

signal_rfile = ROOTFile(datadir("sims","SDBDRC_vertex_prob","RH037_foil_bulk_SDBDRC_vertex_prob.root"))
signal_data = LazyTree(signal_rfile, "tree", ["reconstructedEnergy1", "reconstructedEnergy2"]) |> DataFrame


signal = Process(
    signal_data.reconstructedEnergy1 + signal_data.reconstructedEnergy2,
    "RH_foil_bulk",
    true,    
    1.0,
    9.019296e7,
    1e8,
    0:100:3500,
    6.07,
)
background = [
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
]

set_signal!(background[1], false)

α = 1.64485362695147
t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)




