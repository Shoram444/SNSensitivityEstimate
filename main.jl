using DrWatson, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

Bi214file = ROOTFile("data/sims/Bi214_1e8E_EnePhiDist_J23.root")
Tl208file = ROOTFile("data/sims/Tl208_1e8E_EnePhiDist_J23.root")
K40file = ROOTFile("data/sims/K40_1e8E_EnePhiDist_J23.root")
Pa234mfile = ROOTFile("data/sims/Pa234m_1e8E_EnePhiDist_J23.root")
bbfile = ROOTFile("data/sims/2nubb_1e8E_EnePhiDist_J23.root")   # Standard 2nubb


Bi214 = fill_from_root_file(Bi214file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])  
Tl208 = fill_from_root_file(Tl208file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]) 
K40 = fill_from_root_file(K40file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]) 
Pa234m = fill_from_root_file(Pa234mfile, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]) 
bb = fill_from_root_file(bbfile, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]) 

###########################################
###########################################
#### Sum energy spectra
Bi214SumE    = Process(Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2, "Bi214", :false, BkgActivityParams[:Bi214], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Tl208SumE    = Process(Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2, "Tl208", :false, BkgActivityParams[:Tl208], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
K40SumE      = Process(K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2, "K40", :false, BkgActivityParams[:K40], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Pa234mSumE   = Process(Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2, "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
bbSumE       = Process(bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2, "2nubb", :true, SigActivityParams[:bb2Standard], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])

with(
    gr;
    xlabel = "min sum energy [keV]", 
    ylabel = "max sum energy [keV]", 
    c = :coolwarm, 
    size = (1600, 1000), 
    thickness_scaling = 1.4, 
    right_margin = 16Plots.mm
)   do

    plot(Bi214SumE.efficiency, title ="sum E Efficiency for "* Bi214SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Bi214SumE, allowedtypes = (String, StepRange), "png")), current())

    plot(Tl208SumE.efficiency, title ="sum E Efficiency for "* Tl208SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Tl208SumE, allowedtypes = (String, StepRange), "png")), current())

    plot(K40SumE.efficiency, title ="sum E Efficiency for "* K40SumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", K40SumE, allowedtypes = (String, StepRange), "png")), current())

    plot(Pa234mSumE.efficiency, title ="sum E Efficiency for "* Pa234mSumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", Pa234mSumE, allowedtypes = (String, StepRange), "png")), current())

    plot(bbSumE.efficiency, title ="sum E Efficiency for "* bbSumE.isotopeName)
    safesave(plotsdir("SumE", savename("efficiency_sumE_", bbSumE, allowedtypes = (String, StepRange), "png")), current())
end

stbSum = get_sToBRatio(bbSumE, Bi214SumE, Tl208SumE, Pa234mSumE, K40SumE)
best_stbSum= get_max_bin(stbSum)
plot(
    stbSum;
    c=:coolwarm, 
    title ="signal to background ratio for SumE channel 2nubb 
            max_stb = $(best_stbSum[:maxBinCount] |> round) at ($(best_stbSum[:minBinEdge]) - $(best_stbSum[:maxBinEdge]))keV"
)
safesave(plotsdir("SumE", savename("StoB_SumE_binning",bbSumE ,"png")), current())

get_estimated_bkg_counts(best_stbSum, SNparams, Bi214SumE, Tl208SumE, Pa234mSumE, K40SumE)

###########################################
###########################################
###########################################
###########################################
###########################################
#### Single energy spectra
Bi214SingleE    = Process(vcat(Bi214.reconstructedEnergy1,  Bi214.reconstructedEnergy2), "Bi214", :false, BkgActivityParams[:Bi214], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Tl208SingleE    = Process(vcat(Tl208.reconstructedEnergy1,  Tl208.reconstructedEnergy2), "Tl208", :false, BkgActivityParams[:Tl208], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
K40SingleE      = Process(vcat(K40.reconstructedEnergy1,  K40.reconstructedEnergy2), "K40", :false, BkgActivityParams[:K40], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Pa234mSingleE   = Process(vcat(Pa234m.reconstructedEnergy1,  Pa234m.reconstructedEnergy2), "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
bbSingleE       = Process(vcat(bb.reconstructedEnergy1,  bb.reconstructedEnergy2), "2nubb", :true, SigActivityParams[:bb2Standard], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])

with(
    xlabel = "min single energy [keV]", 
    ylabel = "max single energy [keV]", 
    c = :coolwarm, size = (1600, 1000), 
    thickness_scaling = 1.4, 
    right_margin = 16Plots.mm
)   do
    plot(Bi214SingleE.efficiency, title ="single E Efficiency for "* Bi214SingleE.isotopeName)
    safesave(plotsdir("SingleE", savename("efficiency_singleE_", Bi214SingleE, allowedtypes = (String, StepRange), "png")), current())

    plot(Tl208SingleE.efficiency, title ="single E Efficiency for "* Tl208SingleE.isotopeName)
    safesave(plotsdir("SingleE", savename("efficiency_singleE_", Tl208SingleE, allowedtypes = (String, StepRange), "png")), current())

    plot(K40SingleE.efficiency, title ="single E Efficiency for "* K40SingleE.isotopeName)
    safesave(plotsdir("SingleE", savename("efficiency_singleE_", K40SingleE, allowedtypes = (String, StepRange), "png")), current())

    plot(Pa234mSingleE.efficiency, title ="single E Efficiency for "* Pa234mSingleE.isotopeName)
    safesave(plotsdir("SingleE", savename("efficiency_singleE_", Pa234mSingleE, allowedtypes = (String, StepRange), "png")), current())

    plot(bbSingleE.efficiency, title ="single E Efficiency for "* bbSingleE.isotopeName)
    safesave(plotsdir("SingleE", savename("efficiency_singleE_", bbSingleE, allowedtypes = (String, StepRange), "png")), current())
end

stbSingle = get_sToBRatio(bbSingleE, Bi214SingleE, Tl208SingleE, Pa234mSingleE)
best_stbSingle= get_max_bin(stbSingle)
plot(
    stbSingle; 
    c=:coolwarm, 
    title = "signal to background ratio for SingleE channel 2nubb 
    max_stb = $(best_stbSingle[:maxBinCount] |> round) at ($(best_stbSingle[:minBinEdge]) - $(best_stbSingle[:maxBinEdge]))keV"
)
safesave(plotsdir("SingleE",savename("StoB_singleE_binning",bbSingleE,"png")), current())

###########################################
###########################################
###########################################
###########################################
###########################################
#### Angular Study
Bi214Phi    = Process(Bi214.phi, "Bi214", :false, BkgActivityParams[:Bi214], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
Tl208Phi    = Process(Tl208.phi, "Tl208", :false, BkgActivityParams[:Tl208], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
K40Phi      = Process(K40.phi, "K40", :false, BkgActivityParams[:K40], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
Pa234mPhi   = Process(Pa234m.phi, "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
bbPhi       = Process(bb.phi, "2nubb", :true, SigActivityParams[:bb2Standard], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])



with(
    xlabel = "min angle", 
    ylabel = "max angle", 
    c = :coolwarm, 
    size = (1600, 1000), 
    thickness_scaling = 1.4, 
    right_margin = 16Plots.mm
)   do
    plot(Bi214Phi.efficiency, title ="Angular Efficiency for "* Bi214Phi.isotopeName)
    safesave(plotsdir("Angular", savename("efficiency_anglular_",Bi214Phi, "png")), current())

    plot(Tl208Phi.efficiency, title ="Angular Efficiency for "* Tl208Phi.isotopeName)
    safesave(plotsdir("Angular", savename("efficiency_anglular_",Tl208Phi, "png")), current())

    plot(K40Phi.efficiency, title ="Angular Efficiency for "* K40Phi.isotopeName)
    safesave(plotsdir("Angular", savename("efficiency_anglular_", K40Phi ,"png")), current())

    plot(Pa234mPhi.efficiency, title ="Angular Efficiency for "* Pa234mPhi.isotopeName)
    safesave(plotsdir("Angular", savename("efficiency_anglular_", Pa234mPhi, "png")), current())

    plot(bbPhi.efficiency, title ="Angular Efficiency for "* bbPhi.isotopeName)
    safesave(plotsdir("Angular", savename("efficiency_anglular_", bbPhi, "png")), current())
end

stbPhi = get_sToBRatio(bbPhi, Bi214Phi, Tl208Phi, Pa234mPhi)
best_stbPhi= get_max_bin(stbPhi)
plot(
    stbPhi;
    c=:coolwarm, 
    title ="signal to background ratio for Angular channel 2nubb 
    max_stb = $(best_stbPhi[:maxBinCount] |> round) at ($(best_stbPhi[:minBinEdge]) - $(best_stbPhi[:maxBinEdge]))keV"
)
safesave(plotsdir("Angular", savename("StoB_Angular_",bbPhi,"png")), current())


