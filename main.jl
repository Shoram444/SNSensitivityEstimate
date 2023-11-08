using DrWatson, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
@quickactivate "SensitivityEstimate"

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



###########################################
###########################################
###########################################
###########################################
###########################################
#### Sum energy spectra
sumEParams = Dict(
    :binning => 0:100:3500,
    :nTotalSim => 1e8
)

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x
    :Bi214 => 1.50 / 1000 * SNparams["m"], # [mBq/kg] converted to [Bq]
    :Tl208 => 0.39 / 1000 * SNparams["m"], # [mBq/kg] converted to [Bq]
    :Pa234m => 17.3 / 1000 * SNparams["m"], # [mBq/kg] converted to [Bq]
    :K40 => 58.7 / 1000 * SNparams["m"] # [mBq/kg] converted to [Bq]
)

Bi214SumE    = Process(Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2, "Bi214", :false, BkgActivityParams[:Bi214], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Tl208SumE    = Process(Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2, "Tl208", :false, BkgActivityParams[:Tl208], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
K40SumE      = Process(K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2, "K40", :false, BkgActivityParams[:K40], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Pa234mSumE   = Process(Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2, "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
bbSumE       = Process(bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2, "2nubb", :true, 1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])

let 
    default(xlabel = "min sum energy [keV]", ylabel = "max sum energy [keV]", c = :coolwarm, size = (1600, 1000), thickness_scaling = 1.4, right_margin = 16Plots.mm)
    plot_efficiency(Bi214SumE, sumEParams[:binning], title ="sum E Efficiency for "* Bi214SumE.isotopeName)
    savefig(plotsdir(savename("efficiency_sumE_", Bi214SumE, allowedtypes = (String, StepRange), "png")))
    # $(Bi214SumE.isotopeName)_$(step(sumEParams[:binning]))-degree_binning.png")

    plot_efficiency(Tl208SumE, sumEParams[:binning], title ="sum E Efficiency for "* Tl208SumE.isotopeName)
    savefig(plotsdir(savename("efficiency_sumE_", Tl208SumE, allowedtypes = (String, StepRange), "png")))
    # $(Tl208SumE.isotopeName)_$(step(sumEParams[:binning]))-degree_binning.png")

    plot_efficiency(K40SumE, sumEParams[:binning], title ="sum E Efficiency for "* K40SumE.isotopeName)
    savefig(plotsdir(savename("efficiency_sumE_", K40SumE, allowedtypes = (String, StepRange), "png")))
    # $(K40SumE.isotopeName)_ph$(step(sumEParams[:binning]))-degree_binning.png")

    plot_efficiency(Pa234mSumE, sumEParams[:binning], title ="sum E Efficiency for "* Pa234mSumE.isotopeName)
    savefig(plotsdir(savename("efficiency_sumE_", Pa234mSumE, allowedtypes = (String, StepRange), "png")))
    # $(Pa234mSumE.isotopeName)_phi$(step(sumEParams[:binning]))-degree_binning.ng"))

    plot_efficiency(bbSumE, sumEParams[:binning], title ="sum E Efficiency for "* bbSumE.isotopeName)
    savefig(plotsdir(savename("efficiency_sumE_", bbSumE, allowedtypes = (String, StepRange), "png")))
    # $(bbSumE.isotopeName)_phi$(step(sumEParams[:binning]))-degree_binning.png")
end

stbSum = @show get_sToBRatio(bbSumE, Bi214SumE, Tl208SumE, Pa234mSumE)
best_stbSum= maximum(stbSum) |> round
minSumE, maxSumE = bbSumE.bins[argmax(stbSum)[1]], bbSumE.bins[argmax(stbSum)[2]]
plot_sToBRatio(bbSumE, Bi214SumE, Tl208SumE, Pa234mSumE; c=:coolwarm, title ="signal to background ratio for SumE channel 2nubb \nmax_stb = $best_stbSum at ($minSumE - $maxSumE)keV")
savefig(plotsdir(savename("StoB_SumE_binning",bbSumE.bins,"png")))


###########################################
###########################################
###########################################
###########################################
###########################################
#### Single energy spectra
singleEParams = Dict(
    :binning => 0:100:3500,
    :nTotalSim => 2e8 # because we have 2 electrons per event
)

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x
    :Bi214 => 1.50 / 1000 * SNparams["m"], # [mBq/kg] converted to [Bq]
    :Tl208 => 0.39 / 1000 * SNparams["m"], # [mBq/kg] converted to [Bq]
    :Pa234m => 17.3 / 1000 * SNparams["m"], # [mBq/kg] converted to [Bq]
    :K40 => 58.7 / 1000 * SNparams["m"] # [mBq/kg] converted to [Bq]
)

Bi214SingleE    = Process(vcat(Bi214.reconstructedEnergy1,  Bi214.reconstructedEnergy2), "Bi214", :false, BkgActivityParams[:Bi214], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Tl208SingleE    = Process(vcat(Tl208.reconstructedEnergy1,  Tl208.reconstructedEnergy2), "Tl208", :false, BkgActivityParams[:Tl208], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
K40SingleE      = Process(vcat(K40.reconstructedEnergy1,  K40.reconstructedEnergy2), "K40", :false, BkgActivityParams[:K40], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Pa234mSingleE   = Process(vcat(Pa234m.reconstructedEnergy1,  Pa234m.reconstructedEnergy2), "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
bbSingleE       = Process(vcat(bb.reconstructedEnergy1,  bb.reconstructedEnergy2), "2nubb", :true, 1, SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])

let 
    default(xlabel = "min single energy [keV]", ylabel = "max single energy [keV]", c = :coolwarm, size = (1600, 1000), thickness_scaling = 1.4, right_margin = 16Plots.mm)
    plot_efficiency(Bi214SingleE, singleEParams[:binning], title ="single E Efficiency for "* Bi214SingleE.isotopeName)
    savefig(plotsdir(savename("efficiency_singleE_", Bi214SingleE, allowedtypes = (String, StepRange), "png")))
    # $(Bi214SingleE.isotopeName)_$(step(singleEParams[:binning]))-degree_binning.png")

    plot_efficiency(Tl208SingleE, singleEParams[:binning], title ="single E Efficiency for "* Tl208SingleE.isotopeName)
    savefig(plotsdir(savename("efficiency_singleE_", Tl208SingleE, allowedtypes = (String, StepRange), "png")))
    # $(Tl208SingleE.isotopeName)_$(step(singleEParams[:binning]))-degree_binning.png")

    plot_efficiency(K40SingleE, singleEParams[:binning], title ="single E Efficiency for "* K40SingleE.isotopeName)
    savefig(plotsdir(savename("efficiency_singleE_", K40SingleE, allowedtypes = (String, StepRange), "png")))
    # $(K40SingleE.isotopeName)_ph$(step(singleEParams[:binning]))-degree_binning.png")

    plot_efficiency(Pa234mSingleE, singleEParams[:binning], title ="single E Efficiency for "* Pa234mSingleE.isotopeName)
    savefig(plotsdir(savename("efficiency_singleE_", Pa234mSingleE, allowedtypes = (String, StepRange), "png")))
    # $(Pa234mSingleE.isotopeName)_phi$(step(singleEParams[:binning]))-degree_binning.ng"))

    plot_efficiency(bbSingleE, singleEParams[:binning], title ="single E Efficiency for "* bbSingleE.isotopeName)
    savefig(plotsdir(savename("efficiency_singleE_", bbSingleE, allowedtypes = (String, StepRange), "png")))
    # $(bbSingleE.isotopeName)_phi$(step(singleEParams[:binning]))-degree_binning.png")
end

stbSingle = @show get_sToBRatio(bbSingleE, Bi214SingleE, Tl208SingleE, Pa234mSingleE)
best_stbSingle= maximum(stbSingle) |> round
minSingleE, maxSingleE = bbSingleE.bins[argmax(stbSingle)[1]], bbSingleE.bins[argmax(stbSingle)[2]]
plot_sToBRatio(bbSingleE, Bi214SingleE, Tl208SingleE, Pa234mSingleE; c=:coolwarm, title ="signal to background ratio for 2nubb \nmax_stb = $best_stbSingle at ($minSingleE - $maxSingleE)keV")
savefig(plotsdir(savename("StoB_singleE_binning",bbSingleE.bins,"png")))

###########################################
###########################################
###########################################
###########################################
###########################################
#### Angular Study
phiParams = Dict(
    :binning => 0:5:180,
    :nTotalSim => 1e8
)


Bi214Phi    = Process(Bi214.phi, "Bi214", :false, BkgActivityParams[:Bi214], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
Tl208Phi    = Process(Tl208.phi, "Tl208", :false, BkgActivityParams[:Tl208], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
K40Phi      = Process(K40.phi, "K40", :false, BkgActivityParams[:K40], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
Pa234mPhi   = Process(Pa234m.phi, "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], phiParams[:nTotalSim], phiParams[:binning])
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

stbPhi = @show get_sToBRatio(bbPhi, Bi214Phi, Tl208Phi, Pa234mPhi)
best_stbPhi= maximum(stbPhi) |> round
minPhi, maxPhi = bbPhi.bins[argmax(stbPhi)[1]], bbPhi.bins[argmax(stbPhi)[2]]
plot_sToBRatio(bbPhi, Bi214Phi, Tl208Phi, Pa234mPhi; c=:coolwarm, title ="signal to background ratio for angular for 2nubb \nmax_stb = $best_stbPhi at ($minPhi - $maxPhi)°")
savefig(plotsdir(savename("StoB_Angular_binning",bbPhi.bins,"png")))

