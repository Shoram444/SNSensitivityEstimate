using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using ColorSchemes, SNSensitivityEstimate, CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings,StatsBase, FHist, Distributions

vertex_energy_tof_roi = Dict(
    :r => (0,100), 
    :sumE => (0, 3500),
    :phi => (0, 180),
    :Pint => (0.04, 1),
    :Pext => (0, 0.01),
    :trackLength1 => (0, 150000),
    :trackLength2 => (0, 150000),
)

backgrounds = [
    # 2nubb
    "bb_foil_bulk", 

    # internal contamination
    "Bi214_foil_bulk", 
    "K40_foil_bulk", 
    "Pa234m_foil_bulk", 
    "Tl208_foil_bulk", 

    # radon
    "Bi210_wire_surface",     
    "Bi214_wire_surface", 

    # detector contamination
    "Bi214_PMT_glass_bulk",
    "Bi214_PMT5_glass_bulk",
    "K40_PMT_glass_bulk", 
    "K40_PMT5_glass_bulk", 
    "Tl208_PMT_glass_bulk", 
    "Tl208_PMT5_glass_bulk", 
    "Bi214_calo_8inch_scin_bulk", 
    "K40_calo_8inch_scin_bulk", 

    # gamma
    "Tl208_hall_bulk", 
    "Bi214_hall_bulk", 
    "K40_hall_bulk",

]



analysisDict    = Dict(:mode => "sumE", :cuts => "vertex_energy_tof_roi")

if( analysisDict[:cuts] == "nemo3_roi" )
    roi = nemo3_roi
elseif( analysisDict[:cuts] == "topology_roi" )
    roi = topology_roi
elseif( analysisDict[:cuts] == "tof_roi" )
    roi = tof_roi
elseif( analysisDict[:cuts] == "vertex_tof_roi" )
    roi = vertex_tof_roi
elseif( analysisDict[:cuts] == "vertex_energy_tof_roi" )
    roi = vertex_energy_tof_roi
elseif( analysisDict[:cuts] == "vertex_length_energy_tof_roi" )
    roi = vertex_lenfth_energy_tof_roi
end

binning = (0, 180, 5)
simdir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41")
mode = analysisDict[:mode] 
cat_colors = ["#4E79A7","#C7E6E3","#9AD0CB","#76B7B2","#4E8F8B","#F28E2B","#C66A00","#E15759","#B22222", "#7BC36D", "#59A14F", "#3F7F3A"]


all_processes = load_data_processes(
    simdir,
    String(mode);
    fwhm = 0.0,
    roi  = roi,
)

signal = get_process("bb0nu_foil_bulk", all_processes) |> first
background = [get_process(b, all_processes) |> first for b in backgrounds]

set_signal!(background[1], false)

Q_keV = SNparams["Q"]
α = 1.64485362695147
constants_to_scale_thalf_to_s_to_b = inv(log2(2)*SNparams["Nₐ"]*SNparams["a"]*SNparams["foilMass"]*SNparams["t"]/SNparams["W"])


t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α; approximate="table"), sigdigits=3)

with_theme(theme_latexfonts()) do
    f = Figure(size = (1400, 900),  fontsize = 36, figure_padding = 16)
    a = Axis(
        f[1,1], 
        xlabel = L"$\varphi_l$ (keV) $$", 
        ylabel = L"$\varphi_u$ (keV) $$", 
        title = L"2D sensitivity map for $0\nu\beta\beta$ for $17.5$kgyr exposure, using $\varphi$ channel", 
        )
    p = plot!(a, t12MapESum)
    # text!(a, 2000, 500, text=lbl)
    Colorbar(f[1,2], p, label=L"sensitivity (yr) $$", scale=identity)
    save(joinpath(scriptsdir("phd_final/07_sensitivity/sensitivity_1D/figs"),
        "nu0bb_mode_$(mode)_cuts_$(analysisDict[:cuts])_ThalfbbESum_$(ThalfbbESum)_kgyr.png"), f)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1400, 900),  fontsize = 36, figure_padding = 16)
    a = Axis(
        f[1,1], 
        xlabel = L"$\varphi_l$ (keV) $$", 
        ylabel = L"$\varphi_u$ (keV) $$", 
        title = L"2D efficiency map for $0\nu\beta\beta$, using $\varphi$ channel", 
        )
    p = plot!(a, signal.efficiency)
    # text!(a, 2000, 500, text=lbl)
    Colorbar(f[1,2], p, label=L"$\varepsilon$(ROI)", scale=identity)
    save(joinpath(scriptsdir("phd_final/07_sensitivity/sensitivity_1D/figs"),
        "nu0bb_mode_$(mode)_efficiency_map.png"), f)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1400, 900),  fontsize = 36, figure_padding = 26)
    a = Axis(
        f[1,1], 
        xlabel = L"$\varphi$ (°) $$", 
        ylabel = L"efficiency (%/$5^{\circ}$)", 
        title = L"$\varepsilon(ROI)$ of $0\nu\beta\beta$ in SuperNEMO", 
        limits = (0, 180, 0, nothing),
    )
    h = Hist1D(signal.dataVector; binedges = signal.bins)
    h.bincounts .= h.bincounts .* inv(signal.nTotalSim) .* 100
    plot!(a, h, color = cat_colors[1], label = L"MC reconstructed $0\nu\beta\beta$")

    # mark Q-value at Q_keV
    # vlines!(a, [Q_keV]; color = :red, linestyle = :dash, linewidth = 3,
        # label = L"$Q_{\beta\beta}$ = %$(Q_keV) keV")
    axislegend(a; position = :lt, patchsize = (45,30), fontsize = 24)
    
    # text!(a, 2000, 500, text=lbl)
    save(joinpath(scriptsdir("phd_final/07_sensitivity/sensitivity_1D/figs"),
        "nu0bb_mode_$(mode)_efficiency_1D.png"), f)
    f
end
using Measurements

bkgs = [restrict(get_bkg_counts_1D(b), best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge]) for b in background]
measurements = [
    sum(measurement.(bincounts(b), binerrors(b))) for b in bkgs
]

using PrettyTables
df = DataFrame(
    source = backgrounds,
    counts_in_ROI = round.(measurements, digits=2),
)
@subset! df :counts_in_ROI .> 0.0
@transform! df :percentage_of_total = Measurements.value.(round.(:counts_in_ROI ./ sum(:counts_in_ROI) * 100, digits=2))
push!(df, (
    source = "Total",
    counts_in_ROI = round(sum(df.counts_in_ROI), sigdigits=4),
    percentage_of_total = 100.0,
))

open(scriptsdir("phd_final/07_sensitivity/sensitivity_1D/","nu0bb_mode_$(mode)bkgs_table.tex"), "w") do io
    pretty_table(
        io,
        df,
        backend = Val(:latex),
    )
end

h = Hist1D(signal.dataVector; binedges = signal.bins)
h = h * inv(signal.nTotalSim)
eff_roi_hist = restrict(h, 2700, 3100)
mean_eff_roi = sum(measurement.(bincounts(eff_roi_hist), binerrors(eff_roi_hist)))
