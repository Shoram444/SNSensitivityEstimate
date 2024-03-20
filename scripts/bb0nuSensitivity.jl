using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))
include(scriptsdir("LoadData.jl"))

signalProcessName = "bb0nu_foil_bulk"
sigProcessESum = get_process(signalProcessName, probProcessesESum11)

bb2nuProcess = get_process("bb_foil_bulk", probProcessesESum11)
bb2nuProcess.signal = false
bkgProcessesESum = [
    get_process("Tl208_foil_bulk", probProcessesESum11),
    get_process("Bi214_foil_bulk", probProcessesESum11),
    get_process("Bi214_foil_surface", probProcessesESum11),
    get_process("Bi214_wire_surface", probProcessesESum11),
    bb2nuProcess,
]

α = 1.64


t12MapESum      = get_tHalf_map(SNparams, α, sigProcessESum, bkgProcessesESum...)
best_t12ESum    = get_max_bin(t12MapESum)
expBkgESum      = get_bkg_counts_ROI(best_t12ESum, bkgProcessesESum...)
effbb           = lookup(sigProcessESum, best_t12ESum)
ThalfbbESum     = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
      b  = $(round(expBkgESum, sigdigits = 3)) 
      T12 ≥  $(ThalfbbESum) yr 
      ε = $(round(effbb, sigdigits = 3)*100)%"


with(
    gr;
    c=:coolwarm,
    size=(1600, 1000),
    thickness_scaling=1.7,
    right_margin=16Plots.mm,
    left_margin=16Plots.mm,
    bottom_margin=16Plots.mm,
) do
    plot(
        t12MapESum,
        title="Summed energy channel \nsensitivity as a function of ROI for \n$(sigProcessESum.isotopeName) @ 90% C.L",
        xlabel="min sum energy [keV]",
        ylabel="max sum energy [keV]",
        c=:coolwarm,
        colorbar_title="\n" * L"T_{1/2} ~[yr]",
    )
    annotate!([(2000, 700, lbl)])

    safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_bb=$(bb2nuProcess.isotopeName)_11dE.png"), current())
    safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_bb=$(bb2nuProcess.isotopeName)_11dE.pdf"), current())

    current()
end

t = range(0, 6, 60)

get_tHalf1(t, b) = get_tHalf(
    SNparams["W"],
    SNparams["foilMass"],
    SNparams["Nₐ"],
    t,
    SNparams["a"],
    effbb,
    b / SNparams["tYear"] * t ,
    α;
    approximate="formula"
)

plot(
    t,
    get_tHalf1.(t, expBkgESum),
    xlabel="running time [yr]",
    ylabel="sensitivity " * L"T_{1/2}^{0\nu}" * "[yr]",
    title="sensitivity to " * L"0\nu\beta\beta" * " as a function of running time",
    label="assuming b = $(round(expBkgESum, sigdigits =4))/2.5yr ",
    legendtitle = "ΔE = 11% @ 1MeV",
    legendtitlefontsize = 12
)
plot!(
    t,
    get_tHalf1.(t, 0.5),
    xlabel="running time [yr]",
    ylabel="sensitivity " * L"T_{1/2}^{0\nu}" * "[yr]",
    title="sensitivity to " * L"0\nu\beta\beta" * " as a function of running time",
    label="assuming b = $(round(0.5, sigdigits =4))/2.5yr ",
)
hline!([4.6e24], c=:black, label="cupid sensitivity")
safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_vs_time_11dE.png"), current())


plot(
    t,
    get_tHalf1.(t, 1.5),
    xlabel="running time [yr]",
    ylabel="sensitivity " * L"T_{1/2}^{0\nu}" * "[yr]",
    title="sensitivity to " * L"0\nu\beta\beta" * " as a function of running time",
    label="b = 1.5 / 2.5yr ",
    legendtitle = "ΔE = 11% @ 1MeV",
    legendtitlefontsize = 12,
    # fill = get_tHalf1.(t, 0.5),
    # fa = 0.4,
    # fillstyle= :x
)
plot!(
    t,
    get_tHalf1.(t, 0.5),
    label="b = 0.5 / 2.5yr ",
    c= 2,
)
plot!(
    t,
    get_tHalf1.(t, 0.5),
    label="b ∈ (0.5, 1.5) / 2.5yr ",
    c= :black,
    fill = get_tHalf1.(t, 1.5),
    fa = 0.4,
    fillstyle= :x,
    lw = 0
)
hline!([4.6e24], c=:black, label="cupid sensitivity")
savefig(plotsdir("SumE","T12_band.png"))


hhs = get_bkg_counts_1D.(bkgProcessesESum)
sig1d = get_bkg_counts_1D(sigProcessESum)

gr()
SensitivityModule.stackedhist(
    hhs, 
    yscale =:log10, 
    ylims = (1e-1, 1e4), 
    label = reshape([p.isotopeName for p in bkgProcessesESum], 1, length(hhs)),
    legend = :best,
    xlabel = "sum energy [keV]",
    ylabel = "estimated counts " * L"[keV^{-1}]",
    title = "estimated sum spectrum",
    xlims = (0, 3500),.
    c= [2 3 4 5 1],
)
savefig(plotsdir("SumE", "estimated1D_nu0_spectrum.png"))


estimated_counts_hist1D(probFiles11["Bi214_foil_bulk"], 
                        BkgActivityParams[:Bi214_foil_bulk], 
                        SNparams["foilMass"], 
                        SNparams["t"], 
                        0:150:3500, 
                        SimulationParams[:Bi214_foil_bulk]) |> integral


for p in probProcessesESum
    @show print_isotope_details(p)
end
