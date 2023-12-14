### A Pluto.jl notebook ###
# v0.19.35

using Markdown
using InteractiveUtils

# ╔═╡ af92ba46-999c-11ee-11a4-8554ebc58c0d
using DrWatson

# ╔═╡ fcfd9507-542c-4bab-9b8b-5e8a00b6bb55
@quickactivate "SensitivityEstimate"

# ╔═╡ 8678f310-e581-41ec-92f2-d5b4e1461396
begin
	using StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Measurements
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
end

# ╔═╡ cf56ae15-6048-48fb-ae49-a8c7cd7390eb
begin 
	_palette = ["#00a0f9", "#ba3030", "#22ac74", "#707070", "#9452bd", "#d76055"]
	theme(
		:dao;
		size = (1200, 800),
		legend = :best,
		guidefontsize = 16,
		tickfontsize = 12,
		titlefontsize = 16,
		legendfontsize = 12,
		left_margin = 4Plots.mm,
		right_margin = 8Plots.mm,
		top_margin = 4Plots.mm,
		bottom_margin = 6Plots.mm,
		thickness_scaling = 1.4,
		:linewidth => 3,
		dpi = 200,
		widen = :true,
		:markerstrokewidth => 1,
		:markerstrokecolor => :black,
		:palette => _palette
	)
end

# ╔═╡ 8d36f771-6342-4237-a3fc-b10727fc18e0
md"""
### First we load all the `.root` input simulation files:
1. **Background processes** (all generated on source foil): 
+ Bi214
+ Tl208
+ K40
+ Pa234m
2. **Signal process**:
> *refined spectrum $\xi_{31}, \xi_{51}$* with the following parameters:
- ``\xi_{31} = 0.37``
- ``\xi_{51} = 0.14``
- ``K^{2\nu} = -0.66``
"""

# ╔═╡ 3a279499-3b4c-4636-b9a6-6183868c3ab9
begin
	Bi214file = ROOTFile(datadir("sims", "Bi214_EnePhiDistPos_J23.root"))
	Tl208file = ROOTFile(datadir("sims", "Tl208_EnePhiDistPos_J23.root"))
	K40file   = ROOTFile(datadir("sims", "K40_EnePhiDistPos_J23.root"))
	Pa234mfile= ROOTFile(datadir("sims", "Pa234m_EnePhiDistPos_J23.root"))
	bbfile    = ROOTFile(datadir("sims", "Xi037_EnePhiDistPos_J23.root"))  
end

# ╔═╡ d9102629-c1d6-4a05-aae4-af49d9c9ddbb
md"""
### Now we transform the raw inptu data into `DataFrame` format for easier handling.
Here we also import only three variables: ``E_1``, ``E_2`` and ``\phi``
"""

# ╔═╡ 13073a06-8c2b-402b-96ff-566db2f879c4
begin
	Bi214 = fill_from_root_file(
		Bi214file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]
	)  
	
	Tl208 = fill_from_root_file(
		Tl208file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]
	) 
	
	K40 = fill_from_root_file(
		K40file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]
	) 
	
	Pa234m = fill_from_root_file(
		Pa234mfile, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]
	) 
	
	bb = fill_from_root_file(
		bbfile, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"]
	) 
end

# ╔═╡ ed904b40-d638-44a3-8e01-1bd86dd1a919
md"""
### Next we import all the important parameters for the analysis. 
The import itself has already been done in the 3rd cell by the `include(scriptsdit("Params.jl"))` command. Now we analyse what's inside.

The file `Params.jl` contains all the input constants that will be used in the analysis. These are divided into groups:
* **Detector Parameters** (Dict variable `SNParams`): 
  * Nₐ = $(SNparams["Nₐ"]) ``\rightarrow`` Avogadro's number in [1/mol]
  * W = $(SNparams["W"]) ``\rightarrow`` Se82 molar mass in [kg/mol]
  * a = $(SNparams["a"]) ``\rightarrow`` abundance/foil enrichment; check number
  * m = $(SNparams["m"]) ``\rightarrow`` foil mass in [kg]
  * t = $(SNparams["t"]) ``\rightarrow`` measurement time in [s]
  * tYear = $(SNparams["tYear"]) ``\rightarrow`` measurement time in [y]
  * SeThalf =  $(round(SNparams["SeThalf"], sigdigits = 4)) ``\rightarrow`` Se82 half life **in [s]**, results from NEMO-3

* **Activity Parameters** (Dict variables `BkgActivityParams` and `SigActivityParams`):
  * **Background** (Taken from NEMO-3 paper)
    * :Bi214 => $(BkgActivityParams[:Bi214]) Bq/kg
    * :Tl208 => $(BkgActivityParams[:Tl208]) Bq/kg
    * :Pa234m => $(BkgActivityParams[:Pa234m]) Bq/kg
    * :K40 => $(BkgActivityParams[:K40]) Bq/kg
  * **Signal** 
    * :Xi037 => $(SigActivityParams[:Xi037]) Bq/kg # a mock value is used here, the activity in reality is the free parameter of this analysis

* **Histogram Parameters** which depend on the analysis channel:
  * Sum energy => $(sumEParams[:binning])
  * Single energy => $(singleEParams[:binning])
  * Phi => $(phiParams[:binning])

Lastly, for obtaining efficiencies, the histograms are divided by the `nTotalSim` parameter representing the total number of simulated events. In this example, for Sum energy and phi channels ``nTotalSim = `` $(sumEParams[:nTotalSim]) and for the single energy channel ``nTotalSim = `` $(singleEParams[:nTotalSim]).
"""

# ╔═╡ 86932308-47bf-4ec9-bd53-8cc6b132d32b
md"""
### Before starting the actual analysis, let's look at the raw data:
(for now we look only at the sum energy channel)

Below, the **normalized Sum Energy spectra** for each process are depicted. It can be seen that there is quite a bit of overlap over the whole range. However, the picture changes (slightly) if we scale the spectra according to their activities and measurement time, thus obtaining an **estimate** of the real spectrum, which SuperNEMO will measure. (Considering only the depicted processes...)
"""

# ╔═╡ e1015a40-3307-46b4-8099-1cac469350a3
with(

        legend =:best, 
        bins = sumEParams[:binning], 
        thickness_scaling = 1.4,
        lw = 4,
        right_margin=12Plots.mm,
        size =(1200, 800),
        norm = :true,
        title = "Sum energy normalized spectra",
		
    ) do
	
	stephist(
		bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2,
		label = "refined "*L"2\nu\beta\beta", 
		fill =0, 
		fa = 0.7, 
		fillstyle = :cross
	)
	
    stephist!(
		Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2,label = "Bi214", 
	)
	
    stephist!(
		Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2,label = "Tl208", 
	)
	
    stephist!(
		Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2,label = "Pa234m", 
	)
	
    stephist!(
		K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2,
		label = "K40", 			
		xlabel = L"E_{sum}/"*" keV", 
        ylabel = "normalized countrate / $(step(sumEParams[:binning]))keV", 
	)

	current()
end

# ╔═╡ 86db346d-5e02-444e-aa1a-892e64bc42b3
md"""
### Now the estimated total spectrum as would be seen by SuperNEMO:
The estimated **total** counts for the full spectrum for each process can be calculated as:
 ``$n = A * t * m * \varepsilon$``

With the previously summarized parameters this is:
* Bi214: n = $(BkgActivityParams[:Bi214].val * SNparams["t"] * SNparams["m"] *		(nrow(Bi214) / sumEParams[:nTotalSim]) |> round)
* Tl208: n = $(BkgActivityParams[:Tl208].val * SNparams["t"] * SNparams["m"] *		(nrow(Tl208) / sumEParams[:nTotalSim]) |> round)
* Pa234m: n = $(BkgActivityParams[:Pa234m].val * SNparams["t"] * SNparams["m"] *		(nrow(Pa234m) / sumEParams[:nTotalSim]) |> round)
* K40: n = $(BkgActivityParams[:K40].val * SNparams["t"] * SNparams["m"] *		(nrow(K40) / sumEParams[:nTotalSim]) |> round)
* Xi31: n = $(SigActivityParams[:Xi037].val * SNparams["t"] * SNparams["m"] *		(nrow(bb) / sumEParams[:nTotalSim]) |> round)

Notice that the main contributions to the background for $2\nu\beta\beta$ come from K40 and Pa234m processes mainly in the lower energy region. 
"""

# ╔═╡ cdf251a6-dc5e-4b21-8789-1daafbb676a5
begin 
	# First we normalize histograms to probability
	h1bbSum = normalize(
		StatsBase.fit( Histogram{Float64},
			bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2 ,
			sumEParams[:binning]
		), mode = :probability
	)
	
	h1Bi214Sum = normalize(
			StatsBase.fit( Histogram{Float64},
			Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2,
			sumEParams[:binning]
		), mode = :probability
	)
	
    h1Tl208Sum = normalize(
		StatsBase.fit( Histogram{Float64},
			Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2,
			sumEParams[:binning]
		), mode = :probability
	)
	
    h1Pa234mSum = normalize(
		StatsBase.fit( Histogram{Float64},
			Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2,
			sumEParams[:binning]
		), mode =:probability
	)
	
    h1K40Sum = normalize(
			StatsBase.fit( Histogram{Float64},
			K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2,
			sumEParams[:binning]
		), mode =:probability
	)

	#scale their bincounts to the activity*time*mass to obtain estimated total counts
	
	h1bbSum.weights = h1bbSum.weights .*(
		SigActivityParams[:Xi037].val * SNparams["t"] * SNparams["m"] *
		(nrow(bb) / sumEParams[:nTotalSim])
	)

	# to get stacked histogram, we add the bincontents of the previous histogram to the next
	h1Bi214Sum.weights = h1bbSum.weights .+ h1Bi214Sum.weights .*( 
		BkgActivityParams[:Bi214].val * SNparams["t"] * SNparams["m"] *
		(nrow(Bi214) / sumEParams[:nTotalSim])
	)
	
	h1Tl208Sum.weights = h1Bi214Sum.weights .+  h1Tl208Sum.weights .*(
		BkgActivityParams[:Tl208].val * SNparams["t"] * SNparams["m"] *
		(nrow(Tl208) / sumEParams[:nTotalSim])
	)
	
	h1Pa234mSum.weights = h1Tl208Sum.weights .+ h1Pa234mSum.weights .*(
		BkgActivityParams[:Pa234m].val * SNparams["t"] * SNparams["m"] *
		(nrow(Pa234m) / sumEParams[:nTotalSim])
	)
	
	h1K40Sum.weights = h1Pa234mSum.weights .+  h1K40Sum.weights .*(
		BkgActivityParams[:K40].val * SNparams["t"] * SNparams["m"] *
		(nrow(K40) / sumEParams[:nTotalSim])
	)

	# plot 
	plot(
		h1bbSum, 
		st= :step, f = 0, fa = 0.7, 
		fillstyle =:cross, label = "refined "*L"2\nu\beta\beta"
	)
	plot!(h1Bi214Sum, st =:step, label = "Bi213")
	plot!(h1Tl208Sum, st =:step, label = "Tl208")
	plot!(h1Pa234mSum, st =:step, label = "Pa234m")
	plot!(
		h1K40Sum, st =:step, label = "K40", 
		xlabel =L"E_{sum} "*"keV", 
		ylabel = "estimated counts / $(step(sumEParams[:binning])) keV",
		title = "Stacked histogram:\nEstimated counts total spectrum"
	)
	
end

# ╔═╡ a43d43d8-5801-44ff-ac6d-6c4ca9d538c1
md"""
## Analysis

To perform the analysis as described in `docdb:5833`. We must follow the steps: 
1. calculate the efficiency maps
2. calculate the corresponding s/b ratios
3. determine best ROI
4. extract ``\bar{b}`` and ``\varepsilon`` in the ROI 
5. calculate ``T_{1/2}``

In this analysis, I've created a data-type `Process` which holds the information about the individual studied processes and I've written a number of methods which calculate the efficiency maps, ROIs and T12. The following example shows the pipeline to obtain the T12:
"""

# ╔═╡ 3f78cca6-9419-451e-aedf-40b1db5795e4
md"""
First we instantiate the `Process`'s for each isotope:
"""

# ╔═╡ cb9eec23-3d92-4f52-a0f6-400363db5589
begin
	Bi214SumE = Process(
		Bi214.reconstructedEnergy1 .+ Bi214.reconstructedEnergy2, sumEBi214Params
	)
	
	Tl208SumE = Process(
		Tl208.reconstructedEnergy1 .+ Tl208.reconstructedEnergy2, sumETl208Params
	)
	
	K40SumE = Process(
		K40.reconstructedEnergy1 .+ K40.reconstructedEnergy2, sumEK40Params
	)
	
	Pa234mSumE = Process(
		Pa234m.reconstructedEnergy1 .+ Pa234m.reconstructedEnergy2, sumEPa234mParams
	)
	
	bbSumE = Process(
		bb.reconstructedEnergy1 .+ bb.reconstructedEnergy2, sumEbbParams
	)
	nothing # so there's no cell output
end

# ╔═╡ ab82a10d-15ca-4528-95e2-e5edfa0dd1c5
md"""
### At construction, the efficiency maps are created as well:
"""

# ╔═╡ d9df2537-7a8a-42dd-a689-6f130d97b525
let 
	p1 = plot(bbSumE.efficiency, widen =:false, c =:coolwarm, xrotation=45, title= "Xi31")
	p2 = plot(Bi214SumE.efficiency, widen =:false, c =:coolwarm, xrotation=45, title= "Bi214")
	p3 = plot(Tl208SumE.efficiency, widen =:false, c =:coolwarm, xrotation=45, title= "Tl208")
	p4 = plot(Pa234mSumE.efficiency, widen =:false, c =:coolwarm, xrotation=45, title= "Pa234m")
	p5 = plot(K40SumE.efficiency, widen =:false, c =:coolwarm, xrotation=45, title= "K40")

	l = @layout [_ a{0.6w} _
				b c
				d e]
	
	plot(p1, p2, p3, p4, p5, 
		size = (1400, 1600),layout = l, 
		xlabel ="E_min", ylabel ="E_max", 
		plot_title = "Efficiency maps", 
		bottom_margins = 1Plots.px,
		top_margins = 1Plots.px,
		left_margins = 1Plots.px,
		right_margins = 20Plots.px,
	)
end

# ╔═╡ ed3a655c-eb27-4cc2-bfb2-6541cd0ca4b5
md"""
#### Combining the **estimated** counts into the efficiency maps we can create the s/b maps:
"""

# ╔═╡ ee5d0838-505b-42e3-89e2-5526de144235
stbSum = get_sToBRatio(bbSumE, Bi214SumE, Tl208SumE, Pa234mSumE, K40SumE);

# ╔═╡ 897e4225-1e99-4840-a4fd-e7dd3bf30415
plot(
	stbSum, 
	c =:coolwarm, widen =:false,
	xlabel ="E_min", ylabel ="E_max",
	title= "signal-over-background map"
)

# ╔═╡ a2510d35-d4c4-4f7a-8ef6-6f1638d91dfc
best_stbSum= get_max_bin(stbSum);

# ╔═╡ ae4f2970-2dc9-4157-9a41-9312d7ce8b12
md"""
### Largest s/b value (best ROI) is found with: 
- ROI: ($(best_stbSum[:minBinEdge]) - $(best_stbSum[:maxBinEdge])) keV
- s/b: $(best_stbSum[:maxBinCount] |> round)
"""

# ╔═╡ d7eb81d3-4ce6-422b-aaff-c8d4c74c20a0
expBkgESum = get_estimated_bkg_counts(best_stbSum, SNparams, Bi214SumE, Tl208SumE, Pa234mSumE, K40SumE);


# ╔═╡ 82ec2f54-2044-498f-81e0-9c8c9816d16c
md"""
### The estimated background count in the ROI is: b = $(expBkgESum |> round)
"""

# ╔═╡ df9d6770-a7f0-4e7f-8599-e50703b114c1
effbb = lookup( bbSumE.efficiency, best_stbSum[:minBinEdge], best_stbSum[:maxBinEdge] );

# ╔═╡ 315e3d2e-a8de-4f56-ac15-0931c3dbfc6f
md"""
### The corresponding signal efficiency in ROI is: $effbb
"""

# ╔═╡ da2d004f-8400-4b90-b4d9-181a1fef318a
md"""
### Finally the calculated sensitivity (assuming Gaussian approximation) for CL = 90% is: 
 ``$T_{1/2}^{2\nu} \geq (1.33 \pm 0.01)\cdot 10^{22} \textrm{yr}$``

"""

# ╔═╡ c5ae4098-ef84-4173-b1d0-d64815e82a0a
ThalfbbESum = get_tHalf(SNparams, effbb, expBkgESum, 1.8)


# ╔═╡ Cell order:
# ╠═af92ba46-999c-11ee-11a4-8554ebc58c0d
# ╠═fcfd9507-542c-4bab-9b8b-5e8a00b6bb55
# ╠═8678f310-e581-41ec-92f2-d5b4e1461396
# ╠═cf56ae15-6048-48fb-ae49-a8c7cd7390eb
# ╟─8d36f771-6342-4237-a3fc-b10727fc18e0
# ╠═3a279499-3b4c-4636-b9a6-6183868c3ab9
# ╟─d9102629-c1d6-4a05-aae4-af49d9c9ddbb
# ╠═13073a06-8c2b-402b-96ff-566db2f879c4
# ╟─ed904b40-d638-44a3-8e01-1bd86dd1a919
# ╟─86932308-47bf-4ec9-bd53-8cc6b132d32b
# ╟─e1015a40-3307-46b4-8099-1cac469350a3
# ╟─86db346d-5e02-444e-aa1a-892e64bc42b3
# ╟─cdf251a6-dc5e-4b21-8789-1daafbb676a5
# ╟─a43d43d8-5801-44ff-ac6d-6c4ca9d538c1
# ╟─3f78cca6-9419-451e-aedf-40b1db5795e4
# ╠═cb9eec23-3d92-4f52-a0f6-400363db5589
# ╟─ab82a10d-15ca-4528-95e2-e5edfa0dd1c5
# ╠═d9df2537-7a8a-42dd-a689-6f130d97b525
# ╟─ed3a655c-eb27-4cc2-bfb2-6541cd0ca4b5
# ╠═ee5d0838-505b-42e3-89e2-5526de144235
# ╠═897e4225-1e99-4840-a4fd-e7dd3bf30415
# ╟─ae4f2970-2dc9-4157-9a41-9312d7ce8b12
# ╠═a2510d35-d4c4-4f7a-8ef6-6f1638d91dfc
# ╟─82ec2f54-2044-498f-81e0-9c8c9816d16c
# ╠═d7eb81d3-4ce6-422b-aaff-c8d4c74c20a0
# ╟─315e3d2e-a8de-4f56-ac15-0931c3dbfc6f
# ╠═df9d6770-a7f0-4e7f-8599-e50703b114c1
# ╟─da2d004f-8400-4b90-b4d9-181a1fef318a
# ╠═c5ae4098-ef84-4173-b1d0-d64815e82a0a
