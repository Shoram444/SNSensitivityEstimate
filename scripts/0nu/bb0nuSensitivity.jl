using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Distributions
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))
include(scriptsdir("LoadData.jl"))

signalProcessName = "bb0nu_foil_bulk"
sigProcessESum = get_process(signalProcessName, probProcessesESum)

bb2nuProcess = get_process("bb_foil_bulk", probProcessesESum)
bb2nuProcess.signal = false
bkgProcessesESum = [
    get_process("Tl208_foil_bulk", probProcessesESum),
    get_process("Bi214_foil_bulk", probProcessesESum),
    get_process("Bi214_foil_surface", probProcessesESum),
    get_process("Bi214_wire_surface", probProcessesESum),
    bb2nuProcess,
]

α = 1.64


t12MapESum = get_tHalf_map(SNparams, α, sigProcessESum, bkgProcessesESum...)
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, bkgProcessesESum...)
effbb = lookup(sigProcessESum, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
      b  = $(round(expBkgESum, sigdigits = 3)) 
      T12 ≥  $(ThalfbbESum) yr 
      ε = $(round(effbb, sigdigits = 3)*100)%"


with(
    pythonplot;
    c=:coolwarm,
    size=(1600, 1000),
    thickness_scaling=1.7,
    right_margin=16Plots.mm,
    left_margin=16Plots.mm,
    bottom_margin=16Plots.mm,
    fontfamily="serif",
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
    # safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_bb=$(bb2nuProcess.isotopeName)_11dE.pdf"), current())

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
    b / SNparams["tYear"] * t,
    α;
    approximate="formula"
)

plot(
    t,
    get_tHalf1.(t, expBkgESum),
    xlabel="running time [yr]",
    ylabel="sensitivity " * L"T_{1/2}^{0\nu}" * "[yr]",
    title="sensitivity to " * L"0\nu\beta\beta" * " as a function of running time",
    label="Frequentist; b = $(round(expBkgESum, sigdigits =4))/2.5yr ",
    legendtitle="ΔE = 11% @ 1MeV",
    legendtitlefontsize=12,
    fontfamily="serif",
    thickness_scaling=1.6,
    legend=:bottomright,
    c=1
)

get_tHalf_Bayess(t, muS=2.35) = log(2) * (effbb * SNparams["foilMass"] * SNparams["a"] * SNparams["Nₐ"] * t) / (muS * SNparams["W"])

hline!(
    [get_tHalf_Bayess(2.5, 2.35)],
    xlabel="running time [yr]",
    ylabel="sensitivity " * L"T_{1/2}^{0\nu}" * "[yr]",
    title="sensitivity to " * L"0\nu\beta\beta" * " as a function of running time",
    label="Bayess SN muS = 2.35 ",
    ls=:dash
)
hline!([4.6e24], c=:black, label="cupid sensitivity")
safesave(plotsdir("SumE", "sensitivity_vs_time_7dE_Bayess.png"), current())



hhs = get_bkg_counts_1D.(bkgProcessesESum)
sig1d = get_bkg_counts_1D(sigProcessESum)

pythonplot()
SensitivityModule.stackedhist(
    hhs,
    yscale=:log10,
    ylims=(1e-4, 1e4),
    label=reshape([p.isotopeName for p in bkgProcessesESum], 1, length(hhs)),
    legend=:best,
    xlabel="sum energy [keV]",
    ylabel="estimated counts " * L"[keV^{-1}]",
    title="estimated sum spectrum",
    xlims=(0, 3500),
    c=[2 3 4 5 1],
    lw=0,
    bar_width=150,
    fontfamily="serif"
)
savefig(plotsdir("SumE", "estimated1D_nu0_spectrum.png"))


estimated_counts_hist1D(
    probFiles["Bi214_foil_bulk"],
    BkgActivityParams[:Bi214_foil_bulk],
    SNparams["foilMass"],
    SNparams["t"],
    0:150:3500,
    SimulationParams[:Bi214_foil_bulk]
) |> integral

bkgNames = ["bb_foil_bulk", "Tl208_foil_surface", "Bi214_foil_bulk", "Bi214_wire_surface", "Bi214_foil_surface"]
bkgFiles = Dict(k => probFiles[k] for k in bkgNames if haskey(probFiles, k))

h = Hist1D(; counttype=Float64, binedges=binningDict[:SumE])
for k in keys(bkgFiles)
    if (k == "bb_foil_bulk")
        a = SigActivityParams[Symbol(k)]
    else
        a = BkgActivityParams[Symbol(k)]
    end
    m = sumEParams[Symbol(k)][:amount]
    t = sumEParams[Symbol(k)][:timeMeas]
    h += estimated_counts_hist1D(bkgFiles[k], a, m, t, binningDict[:SumE], SimulationParams[Symbol(k)])
end

minROI, nBins, bStep = 2500, 9, 100
bins = []
bExp = []
for (i, ROI) in enumerate(minROI:bStep:minROI+(bStep*nBins))
    binning = ROI:bStep:ROI+bStep
    push!(bins, binning)

    b = restrict(h, first(binning), last(binning)) |> integral
    push!(bExp, b)
end
expected_bkg_cts_per_ROI = DataFrame(bins=bins, bExp=bExp)

DataFrames.pretty_table(expected_bkg_cts_per_ROI, backend=Val(:markdown))

b_n1 = @. rand(Poisson(expected_bkg_cts_per_ROI[!, :bExp]))
expected_bkg_cts_per_ROI.b_n1 .= b_n1

expected_bkg_cts_per_ROI

with(gr) do
    bar(2500:100:3500, expected_bkg_cts_per_ROI.b_n1)
end

using PrettyTables
pretty_table(expected_bkg_cts_per_ROI, backend=Val(:markdown))
