using DrWatson, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist

push!(LOAD_PATH, srcdir())
using SensitivityModule
Revise.track(SensitivityModule)

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

#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23, # Avogadro's number in 1/mol
    "W" => 0.08192,     # Se82 molar mass in kg/mol
    "a" => 0.99,        # abundance/foil enrichment; check number
    "m" => 6.25,        # foil mass in kg
    "t" => 2.5,         # measurement time in yr
)

#### Plotting of spectra

# generate_raw_plots(Bi214, "Bi214")
# generate_raw_plots(Tl208, "Tl208")
# generate_raw_plots(K40, "K40")
# generate_raw_plots(Pa234m, "Pa234m")
# generate_raw_plots(bb, "2nubb")

#### Sum energy spectra
sumEParams = Dict(
    :binning => 0:100:3500,
    :nTotalSim => 1e8
)


Bi214SumE    = Process(Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2, "Bi214", :false, 0.1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Tl208SumE    = Process(Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2, "Tl208", :false, 0.1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
K40SumE      = Process(K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2, "K40", :false, 0.1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Pa234mSumE   = Process(Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2, "Pa234m", :false, 0.1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
bbSumE       = Process(bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2, "2nubb", :true, 1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])



let 
    default(xlabel = "min sum energy [keV]", ylabel = "max sum energy [keV]", c = :coolwarm, size = (1600, 1000), thickness_scaling = 1.4, right_margin = 16Plots.mm)
    plot_efficiency(Bi214SumE, sumEParams[:binning], title ="sum E Efficiency for "* Bi214SumE.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_sumE_$(Bi214SumE.isotopeName)_$(step(sumEParams[:binning]))-degree_binning.png"))

    plot_efficiency(Tl208SumE, sumEParams[:binning], title ="sum E Efficiency for "* Tl208SumE.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_sumE_$(Tl208SumE.isotopeName)_$(step(sumEParams[:binning]))-degree_binning.png"))

    plot_efficiency(K40SumE, sumEParams[:binning], title ="sum E Efficiency for "* K40SumE.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_sumE_$(K40SumE.isotopeName)_ph$(step(sumEParams[:binning]))-degree_binning.png"))

    plot_efficiency(Pa234mSumE, sumEParams[:binning], title ="sum E Efficiency for "* Pa234mSumE.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_sumE_$(Pa234mSumE.isotopeName)_phi$(step(sumEParams[:binning]))-degree_binning.png"))

    plot_efficiency(bbSumE, sumEParams[:binning], title ="sum E Efficiency for "* bbSumE.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_sumE_$(bbSumE.isotopeName)_phi$(step(sumEParams[:binning]))-degree_binning.png"))
end


plot_sToBRatio(bbSumE, Bi214SumE, Tl208SumE, Pa234mSumE; c=:coolwarm)
maximum(get_sToBRatio(bbSumE, Bi214SumE, Tl208SumE, Pa234mSumE))


lookup_efficiency(bbPhi, 0, 179)
plot_efficiency(bbPhi, bbPhi.bins)

bbPhi.efficiency

#### Angular Study
phiParams = Dict(
    :binning => 0:5:180,
    :nTotalSim => 1e8
)


Bi214Phi    = Process(Bi214.phi, "Bi214", :false, 0.1, SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
Tl208Phi    = Process(Tl208.phi, "Tl208", :false, 0.1, SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
K40Phi      = Process(K40.phi, "K40", :false, 0.1, SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
Pa234mPhi   = Process(Pa234m.phi, "Pa234m", :false, 0.1, SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
bbPhi       = Process(bb.phi, "2nubb", :true, 1, SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])



let 
    default(xlabel = "min angle", ylabel = "max angle", c = :coolwarm, size = (1600, 1000), thickness_scaling = 1.4, right_margin = 16Plots.mm)
    plot_efficiency(Bi214Phi, phiParams[:binning], title ="Angular Efficiency for "* Bi214Phi.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_anglular_$(Bi214Phi.isotopeName)_$(step(phiParams[:binning]))-degree_binning.png"))

    plot_efficiency(Tl208Phi, phiParams[:binning], title ="Angular Efficiency for "* Tl208Phi.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_anglular_$(Tl208Phi.isotopeName)_$(step(phiParams[:binning]))-degree_binning.png"))

    plot_efficiency(K40Phi, phiParams[:binning], title ="Angular Efficiency for "* K40Phi.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_anglular_$(K40Phi.isotopeName)_ph$(step(phiParams[:binning]))-degree_binning.png"))

    plot_efficiency(Pa234mPhi, phiParams[:binning], title ="Angular Efficiency for "* Pa234mPhi.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_anglular_$(Pa234mPhi.isotopeName)_phi$(step(phiParams[:binning]))-degree_binning.png"))

    plot_efficiency(bb, phiParams[:binning], title ="Angular Efficiency for "* bb.isotopeName)
    savefig(joinpath(plotsdir(), "efficiency_anglular_$(Pa234mPhi.isotopeName)_phi$(step(phiParams[:binning]))-degree_binning.png"))
end


plot_sToBRatio(bbPhi, Bi214Phi, Tl208Phi; c=:coolwarm)

lookup_efficiency(bbPhi, 0, 179)
plot_efficiency(bbPhi, bbPhi.bins)

bbPhi.efficiency