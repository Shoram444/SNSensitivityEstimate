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


bbSumVector = bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2 
PaSumVector = Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2 

bbSumE       = Process(bbSumVector, "2nubb", :true, SigActivityParams[:bb2Standard], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])
Pa234mSumE   = Process(PaSumVector, "Pa234m", :false, BkgActivityParams[:Pa234m], SNparams["t"], sumEParams[:nTotalSim], sumEParams[:binning])

h1bb = Hist1D(bbSumVector, 0:200:2800) |> normalize
h1bb = h1bb * bbSumE.activity * SNparams["t"] * SNparams["m"] * (length(bbSumVector) / 1e8)

h1Pa = Hist1D(PaSumVector, 0:200:2800) |> normalize
h1Pa = h1Pa * Pa234mSumE.activity * SNparams["t"] * SNparams["m"] * (length(PaSumVector) / 1e8)

plot(
    h1bb + h1Pa ,
    lw = 3,
    linecolor = 2,
    c = 2,
    fa = 0.2,
    xlabel = L"E_{sum}"*" [keV]",
    ylabel = "no. events / 200keV",
    label = "signal",
    title ="reconstructed energy spectrum",
    thickness_scaling = 1.6,
    legend =:best
)
safesave(plotsdir("Presentation", "sign_estimated_2process.png"), current())

plot!(
    h1Pa,
    lw = 3,
    linecolor = 4,
    c = 4,
    fa = 0.2,
    label = "background",
    thickness_scaling = 1.6,
    legend =:best

)
safesave(plotsdir("Presentation", "sign_and_bkg_estimated_2process.png"), current())

n = length(0:200:2800)
anim = @animate for i = 1:n-1, j = i:n-1
    plot(
        h1bb + h1Pa ,
        lw = 3,
        linecolor = 2,
        c = 2,
        fa = 0.2,
        xlabel = L"E_{sum}"*" [keV]",
        ylabel = "no. events ",
        label = "signal",
        title ="reconstructed energy spectrum",
        thickness_scaling = 1.6,
        legend =:best
    )
    plot!(
        h1Pa,
        lw = 3,
        linecolor = 4,
        c = 4,
        fa = 0.2,
        thickness_scaling = 1.6,
        legend =:best,
        label = "background",
    )
    plot!(
        bincenters(h1bb + h1Pa)[i:j], 
        bincounts(h1bb + h1Pa)[i:j],
        st = :bar,
        normed = :false,
        lw = 3,
        linecolor = 3,
        c = 3,
        fa = 0.5,
        label = "",
        bar_width= 200,
        thickness_scaling = 1.6,
        legend =:best
    )
    eff = round(sum(bincounts(h1bb)[i:j]) / integral(h1bb + h1Pa), sigdigits = 3)

    annotatewithbox!( 
        current(),
        text("$(i*200)-$(j*200+200) keV
        ε = $eff"),
        2300,
        50,
        900,
        10,
    )
    safesave(plotsdir("Presentation", "ROI1D", "ROI=$(i*200)-$(j*200+200).png"), current())
end
gif(anim, plotsdir("Presentation", "ROI1D", "ROI_15bins_2fps_gif.gif"), fps = 2)

