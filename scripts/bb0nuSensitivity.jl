using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))
include(scriptsdir("LoadData.jl"))

signalProcessName = "bb0nu_foil_bulk"
sigProcessESum = get_process(signalProcessName, probProcessesESum)

bb2nuProcess = get_process("Xi037_foil_bulk", probProcessesESum)
bb2nuProcess.signal = false
bkgProcessesESum = [
    get_process("Tl208_foil_bulk", probProcessesESum),
    get_process("Bi214_foil_bulk", probProcessesESum),
    get_process("Bi214_foil_surface", probProcessesESum),
    get_process("Bi214_wire_surface", probProcessesESum),
    bb2nuProcess
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

    safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_bb=$(bb2nuProcess.isotopeName).png"), current())
    safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_bb=$(bb2nuProcess.isotopeName).pdf"), current())

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
safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName)_vs_time.png"), current())