using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

###################################################
###################################################
#### LOADING FILES
###################################################
###################################################

# hack for quick fill from root
ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])

#Bi214
Bi214_foil_bulk = ROOTFile(datadir("sims/Bi214_foil_bulk_1e8E.root")) |> ffrf
Bi214_foil_surface = ROOTFile(datadir("sims/Bi214_foil_surface_1e7E.root")) |> ffrf
Bi214_PMT_bulk = ROOTFile(datadir("sims/Bi214_PMT_bulk_1e7E.root")) |> ffrf
Bi214_wire_bulk = ROOTFile(datadir("sims/Bi214_wire_bulk_1e7E.root")) |> ffrf
Bi214_wire_surface = ROOTFile(datadir("sims/Bi214_wire_surface_1e7E.root")) |> ffrf

#Tl208
Tl208_foil_bulk = ROOTFile(datadir("sims/Tl208_foil_bulk_1e8E.root")) |> ffrf
Tl208_PMT_bulk = ROOTFile(datadir("sims/Tl208_PMT_bulk_1e8E.root")) |> ffrf

#Pa234m
Pa234m_foil_bulk = ROOTFile(datadir("sims/Pa234m_foil_bulk_1e8E.root")) |> ffrf

#K40
K40_foil_bulk = ROOTFile(datadir("sims/K40_foil_bulk_1e8E.root")) |> ffrf
K40_PMT_bulk = ROOTFile(datadir("sims/K40_PMT_bulk_1e8E.root")) |> ffrf

#Xi31
Xi037_foil_bulk = ROOTFile(datadir("sims/Xi037_foil_bulk_1e8E.root")) |> ffrf


###################################################
###################################################
#### INITIATE PROCESS ESUM
###################################################
###################################################

Bi214_foil_bulk_SumE = Process(Bi214_foil_bulk.reconstructedEnergy1 .+ Bi214_foil_bulk.reconstructedEnergy2, sumEBi214_foil_bulk_Params)
Bi214_foil_surface_SumE = Process(Bi214_foil_surface.reconstructedEnergy1 .+ Bi214_foil_surface.reconstructedEnergy2, sumEBi214_foil_surface_Params)
Bi214_PMT_bulk_SumE = Process(Bi214_PMT_bulk.reconstructedEnergy1 .+ Bi214_PMT_bulk.reconstructedEnergy2, sumEBi214_PMT_bulk_Params)
Bi214_wire_bulk_SumE = Process(Bi214_wire_bulk.reconstructedEnergy1 .+ Bi214_wire_bulk.reconstructedEnergy2, sumEBi214_wire_bulk_Params)
Bi214_wire_surface_SumE = Process(Bi214_wire_surface.reconstructedEnergy1 .+ Bi214_wire_surface.reconstructedEnergy2, sumEBi214_wire_surface_Params)


Tl208_foil_bulk_SumE = Process(Tl208_foil_bulk.reconstructedEnergy1 .+ Tl208_foil_bulk.reconstructedEnergy2, sumETl208_foil_bulk_Params)
Tl208_PMT_bulk_SumE = Process(Tl208_PMT_bulk.reconstructedEnergy1 .+ Tl208_PMT_bulk.reconstructedEnergy2, sumETl208_PMT_bulk_Params)


Pa234m_foil_bulk_SumE = Process(Pa234m_foil_bulk.reconstructedEnergy1 .+ Pa234m_foil_bulk.reconstructedEnergy2, sumEPa234m_foil_bulk_Params)


K40_foil_bulk_SumE = Process(K40_foil_bulk.reconstructedEnergy1 .+ K40_foil_bulk.reconstructedEnergy2, sumEK40_foil_bulk_Params)
K40_PMT_bulk_SumE = Process(K40_PMT_bulk.reconstructedEnergy1 .+ K40_PMT_bulk.reconstructedEnergy2, sumEK40_PMT_bulk_Params)


Xi037_foil_bulk_SumE = Process(Xi037_foil_bulk.reconstructedEnergy1 .+ Xi037_foil_bulk.reconstructedEnergy2, sumEXiParams)

###################################################
###################################################
#### PLOT AND SAVE EFFICIENCY MAPS
###################################################
###################################################
with(
    gr;
    xlabel="min sum energy [keV]",
    ylabel="max sum energy [keV]",
    c=:coolwarm,
    size=(1600, 1000),
    thickness_scaling=1.4,
    right_margin=16Plots.mm,
    left_margin=16Plots.mm,
    bottom_margin=16Plots.mm,
) do
    #Bi214
    plot(Bi214_foil_bulk_SumE.efficiency, title="sum E Efficiency for " * Bi214_foil_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Bi214_foil_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    plot(Bi214_foil_surface_SumE.efficiency, title="sum E Efficiency for " * Bi214_foil_surface_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Bi214_foil_surface_SumE, allowedtypes=(String, StepRange), "png")), current())

    plot(Bi214_wire_bulk_SumE.efficiency, title="sum E Efficiency for " * Bi214_wire_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Bi214_wire_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    plot(Bi214_wire_surface_SumE.efficiency, title="sum E Efficiency for " * Bi214_wire_surface_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Bi214_wire_surface_SumE, allowedtypes=(String, StepRange), "png")), current())

    plot(Bi214_PMT_bulk_SumE.efficiency, title="sum E Efficiency for " * Bi214_PMT_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Bi214_PMT_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())


    #Tl208   
    plot(Tl208_foil_bulk_SumE.efficiency, title="sum E Efficiency for " * Tl208_foil_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Tl208_foil_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    plot(Tl208_PMT_bulk_SumE.efficiency, title="sum E Efficiency for " * Tl208_PMT_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Tl208_PMT_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    #Pa234m
    plot(Pa234m_foil_bulk_SumE.efficiency, title="sum E Efficiency for " * Pa234m_foil_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Pa234m_foil_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    #K40
    plot(K40_foil_bulk_SumE.efficiency, title="sum E Efficiency for " * K40_foil_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", K40_foil_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    plot(K40_PMT_bulk_SumE.efficiency, title="sum E Efficiency for " * K40_PMT_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", K40_PMT_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())

    #Xi31
    plot(Xi037_foil_bulk_SumE.efficiency, title="sum E Efficiency for " * Xi037_foil_bulk_SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Xi037_foil_bulk_SumE, allowedtypes=(String, StepRange), "png")), current())
end


stbSum = get_sToBRatio(
    Bi214_foil_bulk_SumE,
    Bi214_foil_surface_SumE,
    Bi214_PMT_bulk_SumE,
    Bi214_wire_bulk_SumE,
    Bi214_wire_surface_SumE,
    Tl208_foil_bulk_SumE,
    Tl208_PMT_bulk_SumE,
    Pa234m_foil_bulk_SumE,
    K40_foil_bulk_SumE,
    K40_PMT_bulk_SumE,
    Xi037_foil_bulk_SumE,
)
best_stbSum = get_max_bin(stbSum)
expBkgESum = get_estimated_bkg_counts(
    best_stbSum,
    SNparams,
    Bi214_foil_bulk_SumE,
    Bi214_foil_surface_SumE,
    Bi214_PMT_bulk_SumE,
    Bi214_wire_bulk_SumE,
    Bi214_wire_surface_SumE,
    Tl208_foil_bulk_SumE,
    Tl208_PMT_bulk_SumE,
    Pa234m_foil_bulk_SumE,
    K40_foil_bulk_SumE,
    K40_PMT_bulk_SumE,
)
effbb = lookup(Xi037_foil_bulk_SumE.efficiency, best_stbSum[:minBinEdge], best_stbSum[:maxBinEdge])
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, 1.8), sigdigits=3)

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
        stbSum, 
        title= "signal-to-background ratio: " * L"\frac{\varepsilon(ROI)}{\bar{b}(ROI)}", 
        xlabel="min sum energy [keV]",
        ylabel="max sum energy [keV]",
        c = :coolwarm,
        label ="$(best_stbSum[:minBinEdge]) - $(best_stbSum[:maxBinEdge]) keV \ns/b = $(best_stbSum[:maxBinCount] |> round) \nb   = $(expBkgESum |> round) \nT12 ≥  $(ThalfbbESum) yr",
        legend=:bottomright
    )

    safesave(plotsdir("SumE", "signal_to_background.png"), current())
    safesave(plotsdir("SumE", "signal_to_background.pdf"), current())

    current()
end

