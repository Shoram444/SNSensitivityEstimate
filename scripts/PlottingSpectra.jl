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

### Plotting raw spectra!
### Single and Angular unnormed
generate_raw_plots(Bi214, "Bi214", c = 1)
generate_raw_plots(Tl208, "Tl208", c = 2)
generate_raw_plots(K40, "K40", c =3)
generate_raw_plots(Pa234m, "Pa234m", c = 4)
generate_raw_plots(bb, "2nubb", c = 5)

### Single energy
with(
        xlabel = "E_single", 
        ylabel = "normalized countrate / 100keV", 
        legend =:best, 
        bins = singleEParams[:binning], 
        thickness_scaling = 1.4,
        lw = 4,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Single energy normalize spectra"
    ) do
    stephist(vcat(Bi214.reconstructedEnergy1, Bi214.reconstructedEnergy2),label = "Bi214", fill =0, fa = 0.)
    stephist!(vcat(Tl208.reconstructedEnergy1, Tl208.reconstructedEnergy2),label = "Tl208", fill =0, fa = 0.)
    stephist!(vcat(K40.reconstructedEnergy1, K40.reconstructedEnergy2),label = "K40", fill =0, fa = 0.)
    stephist!(vcat(Pa234m.reconstructedEnergy1, Pa234m.reconstructedEnergy2),label = "Pa234m", fill =0, fa = 0.)
    stephist!(vcat(bb.reconstructedEnergy1, bb.reconstructedEnergy2),label = L"2\nu\beta\beta", fill =0, fa = 0.)
    safesave(plotsdir("Raw", "Single", "Normed_Together_single_energy.png"), current())
end

with(
        xlabel = "E_single", 
        ylabel = "Estimated counts / 100keV", 
        legend =:best, 
        bins = singleEParams[:binning], 
        thickness_scaling = 1.4,
        lw = 2,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Estimated spectra based on Activity * Tmeas (2.5yr)"
    ) do
    h1Bi214Single = Hist1D(vcat(Bi214.reconstructedEnergy1, Bi214.reconstructedEnergy1 ),singleEParams[:binning])
    h1Tl208Single = Hist1D(vcat(Tl208.reconstructedEnergy1, Tl208.reconstructedEnergy1 ),singleEParams[:binning])
    h1Pa234mSingle = Hist1D(vcat(Pa234m.reconstructedEnergy1, Pa234m.reconstructedEnergy1 ),singleEParams[:binning])
    h1K40Single = Hist1D(vcat(K40.reconstructedEnergy1, K40.reconstructedEnergy1 ),singleEParams[:binning])
    h1bbSingle = Hist1D(vcat(bb.reconstructedEnergy1, bb.reconstructedEnergy1 ),singleEParams[:binning])
    
    stackHistSingle = stackedhist(
        [h1Bi214Single, h1Tl208Single, h1Pa234mSingle, h1K40Single, h1bbSingle] ,
        [BkgActivityParams[:Bi214], BkgActivityParams[:Tl208], BkgActivityParams[:Pa234m], BkgActivityParams[:K40], SigActivityParams[:bb2Standard]] .* SNparams["t"] .* SNparams["m"],
        label= ["Bi215" "Tl208" "Pa234m" "K40" L"2\nu\beta\beta"],
        # yscale =:log10
    )
    safesave(plotsdir("Raw", "Single", "Estimated_counts_together.png"), current())
    current()
end


### Sum energy
with(
        xlabel = "E_sum", 
        ylabel = "normalized countrate / 100keV", 
        legend =:best, 
        bins = sumEParams[:binning], 
        thickness_scaling = 1.4,
        lw = 2,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Sum energy normalized spectra"
    ) do
    stephist(Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2,label = "Bi214", fill =0, fa = 0.2)
    stephist!(Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2,label = "Tl208", fill =0, fa = 0.2)
    stephist!(K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2,label = "K40", fill =0, fa = 0.2)
    stephist!(Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2,label = "Pa234m", fill =0, fa = 0.2)
    stephist!(bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2,label = L"2\nu\beta\beta", fill =0, fa = 0.2)
    safesave(plotsdir("Presentation", "Sum", "Normed_Together_sum_energy.png"), current())
end

with(
        xlabel = "E_sum", 
        ylabel = "Estimated counts / 100keV", 
        legend =:best, 
        bins = sumEParams[:binning], 
        thickness_scaling = 1.4,
        lw = 3,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Estimated spectra based on Activity * Tmeas (2.5yr)"
    ) do
    h1Bi214Sum = Hist1D(Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy1 ,sumEParams[:binning])
    h1Tl208Sum = Hist1D(Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy1 ,sumEParams[:binning])
    h1Pa234mSum = Hist1D(Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy1 ,sumEParams[:binning])
    h1K40Sum = Hist1D(K40.reconstructedEnergy1 .+ K40.reconstructedEnergy1 ,sumEParams[:binning])
    h1bbSum = Hist1D(bb.reconstructedEnergy1 .+ bb.reconstructedEnergy1 ,sumEParams[:binning])
    
    stackHistSum = stackedhist(
        [h1Bi214Sum, h1Tl208Sum, h1Pa234mSum, h1K40Sum, h1bbSum] ,
        [BkgActivityParams[:Bi214], BkgActivityParams[:Tl208], BkgActivityParams[:Pa234m], BkgActivityParams[:K40], SigActivityParams[:bb2Standard]] .* SNparams["t"] .* SNparams["m"],
        label= ["Bi215" "Tl208" "Pa234m" "K40" L"2\nu\beta\beta"],
        # yscale =:log10
    )
    safesave(plotsdir("Presentation", "Sum", "Estimated_counts_together.png"), current())
    current()
end


### Angular
with(
        xlabel = "angle", 
        ylabel = "normalized countrate / 5°", 
        legend =:best, 
        bins = phiParams[:binning], 
        thickness_scaling = 1.4,
        lw = 2,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Angular normalized spectra"
    ) do
    stephist(Bi214.phi,label = "Bi214", fill =0, fa = 0.2)
    stephist!(Tl208.phi,label = "Tl208", fill =0, fa = 0.2)
    stephist!(K40.phi,label = "K40", fill =0, fa = 0.2)
    stephist!(Pa234m.phi,label = "Pa234m", fill =0, fa = 0.2)
    stephist!(bb.phi,label = L"2\nu\beta\beta", fill =0, fa = 0.2)
    safesave(plotsdir("Raw", "Angular", "Normed_Together_angular.png"), current())
end

with(
        xlabel = "angle", 
        ylabel = "Estimated counts / 5°", 
        legend =:best, 
        bins = phiParams[:binning], 
        thickness_scaling = 1.4,
        lw = 3,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Estimated spectra based on Activity * Tmeas (2.5yr)"
    ) do
    h1Bi214Phi = Hist1D(Bi214.phi ,phiParams[:binning])
    h1Tl208Phi = Hist1D(Tl208.phi ,phiParams[:binning])
    h1Pa234mPhi = Hist1D(Pa234m.phi ,phiParams[:binning])
    h1K40Phi = Hist1D(K40.phi ,phiParams[:binning])
    h1bbPhi = Hist1D(bb.phi ,phiParams[:binning])
    
    stackHistPhi = stackedhist(
        [h1Bi214Phi, h1Tl208Phi, h1Pa234mPhi, h1K40Phi, h1bbPhi] ,
        [BkgActivityParams[:Bi214], BkgActivityParams[:Tl208], BkgActivityParams[:Pa234m], BkgActivityParams[:K40], SigActivityParams[:bb2Standard]] .* SNparams["t"] .* SNparams["m"],
        label= ["Bi215" "Tl208" "Pa234m" "K40" L"2\nu\beta\beta"],
        # yscale =:log10
    )
    safesave(plotsdir("Raw", "Angular", "Estimated_counts_together.png"), current())
    current()
end