using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))
include(scriptsdir("LoadData.jl"))

signalProcessName = "bb0nu_foil_bulk"

###################################################
###################################################
#### INITIATE PROCESS ESUM
###################################################
###################################################

α = 1.64

bkgProcessesESum = [
    p for p in probProcessesESum if !p.signal
]

sigProcessESum = get_process(signalProcessName , probProcessesESum )

t12MapESum = get_tHalf_map(SNparams, α, sigProcessESum, bkgProcessesESum...)
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(
    best_t12ESum,
    bkgProcessesESum...
)
effbb = lookup(sigProcessESum, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)
lbl="$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV \nb  = $(expBkgESum |> round) \nT12 ≥  $(ThalfbbESum) yr"

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
        colorbar_title = "\n" * L"T_{1/2} ~[yr]",
    )
    annotate!([( 2000,700, lbl )])

    safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName).png"), current())
    safesave(plotsdir("SumE", "T12Map_$(sigProcessESum.isotopeName).pdf"), current())

    current()
end

######### SINGLE E

bkgProcessesESingle = [
    p for p in probProcessesESingle if !p.signal
]

sigProcessESingle = get_process(signalProcessName , probProcessesESingle )

t12MapESingle = get_tHalf_map(SNparams, α, sigProcessESingle, bkgProcessesESingle...)
best_t12ESingle = get_max_bin(t12MapESingle)
expBkgESingle = get_bkg_counts_ROI(
    best_t12ESingle,
    bkgProcessesESingle...
)
effbb = lookup(sigProcessESingle, best_t12ESingle)
ThalfbbESingle = round(get_tHalf(SNparams, effbb, expBkgESingle, α), sigdigits=3)
lbl="\n $(best_t12ESingle[:minBinEdge]) - $(best_t12ESingle[:maxBinEdge]) keV \nb   = $(expBkgESingle |> round) \nT12 ≥  $(ThalfbbESingle) yr"

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
        t12MapESingle,
        title="Single electron channel \nsensitivity as a function of ROI for \n$(sigProcessESum.isotopeName) @ 90% C.L",
        xlabel="min single energy [keV]",
        ylabel="max single energy [keV]",
        c=:coolwarm,
        colorbar_title = "\n" * L"T_{1/2} ~[yr]",
    )
    annotate!([( 2000,700, lbl )])

    safesave(plotsdir("SingleE", "T12Map_$(sigProcessESum.isotopeName).png"), current())
    safesave(plotsdir("SingleE", "T12Map_$(sigProcessESum.isotopeName).pdf"), current())

    current()
end


####### phi
bkgProcessesPhi = [
    p for p in probProcessesPhi if !p.signal
]

sigProcessPhi = get_process(signalProcessName , probProcessesPhi )

t12MapPhi = get_tHalf_map(SNparams, α, sigProcessPhi, bkgProcessesPhi...)
best_t12Phi = get_max_bin(t12MapPhi)
expBkgPhi = get_bkg_counts_ROI(
    best_t12Phi,
    bkgProcessesPhi...
)
effbb = lookup(sigProcessPhi, best_t12Phi)
ThalfbbPhi = round(get_tHalf(SNparams, effbb, expBkgPhi, α), sigdigits=3)
lbl="\n $(best_t12Phi[:minBinEdge]) - $(best_t12Phi[:maxBinEdge]) keV \nb   = $(expBkgPhi |> round) \nT12 ≥  $(ThalfbbESingle) yr"

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
        t12MapPhi,
        title="Angular distribution channel \nsensitivity as a function of ROI for \n$(sigProcessESum.isotopeName) @ 90% C.L",
        xlabel="min angle [°]",
        ylabel="max angle [°]",
        c=:coolwarm,
        colorbar_title = "\n" * L"T_{1/2} ~[yr]",
    )
    annotate!([( 120,60, lbl )])

    safesave(plotsdir("Angular", "T12Map_$(sigProcessESum.isotopeName).png"), current())
    safesave(plotsdir("Angular", "T12Map_$(sigProcessESum.isotopeName).pdf"), current())

    current()
end

