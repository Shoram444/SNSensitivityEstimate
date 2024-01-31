### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ af92ba46-999c-11ee-11a4-8554ebc58c0d
using DrWatson

# ╔═╡ fcfd9507-542c-4bab-9b8b-5e8a00b6bb55
@quickactivate "SensitivityEstimate"

# ╔═╡ 8678f310-e581-41ec-92f2-d5b4e1461396
begin
	using StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Measurements, PrettyTables
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
end

# ╔═╡ cf56ae15-6048-48fb-ae49-a8c7cd7390eb
begin 
	_palette = ["#00a0f9", "#ba3030", "#22ac74", "#707070", "#9452bd", "#80ff00", "#ffcc00", "#ff00ff", "#00ffff", "#cc9900", "#2a2727"]
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
		widen = :false,
		:markerstrokewidth => 1,
		:markerstrokecolor => :black,
		:palette => _palette
	)
end

# ╔═╡ 8d36f771-6342-4237-a3fc-b10727fc18e0
md"""
### First we load all the `.root` input simulation files:
1. **Background processes** : 
+ Bi214 foil bulk SumE
+ Bi214 foil surface SumE
+ Bi214 PMT bulk SumE
+ Bi214 wire bulk SumE
+ Bi214 wire surface SumE
+ Tl208 foil bulk SumE
+ Tl208 PMT bulk SumE
+ Pa234m foil bulk SumE
+ K40 foil bulk SumE
+ K40 PMT bulk SumE
2. **Signal process**:
> *refined spectrum $\xi_{31}, \xi_{51}$* with the following parameters:
- ``\xi_{31} = 0.37``
- ``\xi_{51} = 0.14``
- ``K^{2\nu} = -0.66``
"""

# ╔═╡ f93ae5c5-0702-4ec2-9294-7d506153c81e
# hack for quick fill from root
ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"] ) 

# ╔═╡ 3a279499-3b4c-4636-b9a6-6183868c3ab9
begin
#Bi214
	Bi214_foil_bulk = 
		ROOTFile(datadir("sims/Bi214_foil_bulk_1e8E.root")) |> ffrf
	Bi214_foil_surface = 
		ROOTFile(datadir("sims/Bi214_foil_surface_1e7E.root")) |> ffrf
	Bi214_PMT_bulk = 
		ROOTFile(datadir("sims/Bi214_PMT_bulk_1e7E.root")) |> ffrf
	Bi214_wire_bulk = 
		ROOTFile(datadir("sims/Bi214_wire_bulk_1e7E.root")) |> ffrf
	Bi214_wire_surface = 
		ROOTFile(datadir("sims/Bi214_wire_surface_1e7E.root")) |> ffrf

#Tl208
	Tl208_foil_bulk = 
		ROOTFile(datadir("sims/Tl208_foil_bulk_1e8E.root")) |> ffrf
	Tl208_PMT_bulk = 
		ROOTFile(datadir("sims/Tl208_PMT_bulk_1e8E.root")) |> ffrf

#Pa234m
	Pa234m_foil_bulk = 
		ROOTFile(datadir("sims/Pa234m_foil_bulk_1e8E.root")) |> ffrf

#K40
	K40_foil_bulk = 
		ROOTFile(datadir("sims/K40_foil_bulk_1e8E.root")) |> ffrf
	K40_PMT_bulk = 
		ROOTFile(datadir("sims/K40_PMT_bulk_1e8E.root")) |> ffrf

#Xi31
	Xi037_foil_bulk = 
		ROOTFile(datadir("sims/Xi037_foil_bulk_1e8E.root")) |> ffrf 
end

# ╔═╡ d9102629-c1d6-4a05-aae4-af49d9c9ddbb
md"""
### Now we transform the raw inptu data into `DataFrame` format for easier handling.
Here we also import only three variables: ``E_1``, ``E_2`` and ``\phi``
"""

# ╔═╡ ed904b40-d638-44a3-8e01-1bd86dd1a919
md"""
### Next we import all the important parameters for the analysis. 
The import itself has already been done in the 3rd cell by the `include(scriptsdit("Params.jl"))` command. Now we analyse what's inside.

The file `Params.jl` contains all the input constants that will be used in the analysis. These are divided into groups:
* **Detector Parameters** (Dict variable `SNParams`): 
  * Nₐ = $(SNparams["Nₐ"]) ``\rightarrow`` Avogadro's number in [1/mol]
  * W = $(SNparams["W"]) ``\rightarrow`` Se82 molar mass in [kg/mol]
  * a = $(SNparams["a"]) ``\rightarrow`` abundance/foil enrichment; check number
  * foilMass = $(SNparams["foilMass"]) ``\rightarrow`` foil mass in [kg]
  * PMTGlassMass => $(SNparams["PMTGlassMass"]) ``\rightarrow`` PMT glass mass in [kg]
  * wireBulkMass = $(SNparams["wireBulkMass"]) ``\rightarrow`` DUMMY VALUE!! mass of the tracker wires in [kg]
  * t = $(SNparams["t"]) ``\rightarrow`` measurement time in [s]
  * tYear = $(SNparams["tYear"]) ``\rightarrow`` measurement time in [y]
  * SeThalf =  $(round(SNparams["SeThalf"], sigdigits = 4)) ``\rightarrow`` Se82 half life **in [s]**, results from NEMO-3

* **Activity Parameters** (Dict variables `BkgActivityParams` and `SigActivityParams`):
  * **Background** (Taken from NEMO-3 paper, SN measurements, presentations, mock)
    * :Bi214 foil bulk => $(BkgActivityParams[:Bi214_foil_bulk]) Bq/kg
    * :Tl208 foil bulk => $(BkgActivityParams[:Tl208_foil_bulk]) Bq/kg
    * :Pa234m foil bulk => $(BkgActivityParams[:Pa234m_foil_bulk]) Bq/kg
    * :K40 foil bulk => $(BkgActivityParams[:K40_foil_bulk]) Bq/kg
    * :Bi214 PMT bulk => $(round(BkgActivityParams[:Bi214_PMT_bulk], sigdigits=3)) Bq/kg
    * :Tl208 PMT bulk => $(round(BkgActivityParams[:Tl208_PMT_bulk], sigdigits =3)) Bq/kg
    * :K40 PMT bulk => $(round(BkgActivityParams[:K40_PMT_bulk], sigdigits = 3)) Bq/kg
    * :Bi214 wires => $(BkgActivityParams[:Bi214_wire_surface]) Bq/kg

  * **Signal** 
    * :Xi037_foil_bulk => $(SigActivityParams[:Xi037_foil_bulk]) Bq/kg # a mock value is used here, the activity in reality is the free parameter of this analysis

* **Histogram Parameters** which depend on the analysis channel:
  * Sum energy => $(sumEParams[:binning])
  * Single energy => $(singleEParams[:binning])
  * Phi => $(phiParams[:binning])

Lastly, for obtaining efficiencies, the histograms are divided by the `nTotalSim` parameter representing the total number of simulated events. Each process has it's own nTotalSim value associated with it.
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
		Xi037_foil_bulk.reconstructedEnergy1 .+ Xi037_foil_bulk.reconstructedEnergy2,
		label = "refined "*L"2\nu\beta\beta", 
		fill =0, 
		fa = 0.7, 
		fillstyle = :cross
	)
	
    stephist!(
		Bi214_foil_bulk.reconstructedEnergy1 .+ Bi214_foil_bulk.reconstructedEnergy2,label = "Bi214_foil_bulk", 
	)
	stephist!(
		Bi214_foil_surface.reconstructedEnergy1 .+ Bi214_foil_surface.reconstructedEnergy2,label = "Bi214_foil_surface", 
	)
	stephist!(
		Bi214_PMT_bulk.reconstructedEnergy1 .+ Bi214_PMT_bulk.reconstructedEnergy2,label = "Bi214_PMT_bulk", 
	)
	stephist!(
		Bi214_wire_bulk.reconstructedEnergy1 .+ Bi214_wire_bulk.reconstructedEnergy2,label = "Bi214_wire_bulk", 
	)
	stephist!(
		Bi214_wire_surface.reconstructedEnergy1 .+ Bi214_wire_surface.reconstructedEnergy2,label = "Bi214_wire_surface", 
	)
	stephist!(
		Tl208_foil_bulk.reconstructedEnergy1 .+ Tl208_foil_bulk.reconstructedEnergy2,label = "Tl208_foil_bulk", 
	)
	stephist!(
		Tl208_PMT_bulk.reconstructedEnergy1 .+ Tl208_PMT_bulk.reconstructedEnergy2,label = "Tl208_PMT_bulk", 
	)
	stephist!(
		Pa234m_foil_bulk.reconstructedEnergy1 .+ Pa234m_foil_bulk.reconstructedEnergy2,label = "Pa234m_foil_bulk", 
	)
	stephist!(
		K40_PMT_bulk.reconstructedEnergy1 .+ K40_PMT_bulk.reconstructedEnergy2,label = "K40_PMT_bulk", 
	)
	stephist!(
		K40_foil_bulk.reconstructedEnergy1 .+ K40_foil_bulk.reconstructedEnergy2,label = "K40_foil_bulk", 
		xlabel = L"E_{sum} "*"[keV]", 
        ylabel = "normalized countrate " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]", 
		xlims = (0,4500)
		
	)

	safesave(plotsdir("SumE", "Normalized_hist_total_spectrum.pdf"), current())
	safesave(plotsdir("SumE", "Normalized_hist_total_spectrum.png"), current())
	
	current()
end

# ╔═╡ 86db346d-5e02-444e-aa1a-892e64bc42b3
md"""
### Now the estimated total spectrum as would be seen by SuperNEMO:
The estimated **total** counts for the full spectrum for each process can be calculated as:
 ``$n = A * t * m * \varepsilon$``

With the previously summarized parameters this is:
"""

# ╔═╡ 01fdc8b3-e467-401e-8ccb-8397b72d057a
function get_isotope_details( activityParams, SNParams, simParams, isotope, dfData )
	a = activityParams[isotope]
	m = if( occursin("foil_bulk", string(isotope)))
		SNParams["foilMass"]
	elseif( occursin("surface", string(isotope)))
		SNParams["gasVolume"]
	elseif( occursin("PMT_bulk", string(isotope)))
		SNParams["PMTGlassMass"]
	elseif( occursin("wire_bulk", string(isotope)))
		SNParams["wireBulkMass"]
	end
	t = SNParams["t"]
	nTotSim = simParams[isotope]
	ε = nrow(dfData) / nTotSim

	nExpTot = a*m*t*ε

	return a, m, t, ε, nExpTot, nTotSim 
end
	
		

# ╔═╡ a5f981b8-aa70-4d1d-a82b-491feb1a1017
function print_isotope_details( activityParams, SNParams, simParams, isotope, dfData )
	(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( activityParams, SNParams, simParams, isotope, dfData ) 

	a = round(a, sigdigits=3)
	ε = round(ε*100, sigdigits=3)
	m = round(m, sigdigits=3)
	nExpTot = round(nExpTot, sigdigits=3)
	
	"|$isotope | $nExpTot | $ε | $a | $m |"
end

# ╔═╡ 89852571-c350-4f50-84e0-5eeba29b5da7
begin
	isotope_df = DataFrame()
	bkgIsotopes = [
		:Bi214_foil_bulk,
		:Bi214_foil_surface,
		:Bi214_PMT_bulk,
		:Bi214_wire_bulk,
		:Bi214_wire_surface,
		:Tl208_foil_bulk,
		:Tl208_PMT_bulk,
		:Pa234m_foil_bulk,
		:K40_foil_bulk,
		:K40_PMT_bulk,
	]
	bkgDfs = [		
		Bi214_foil_bulk,
		Bi214_foil_surface,
		Bi214_PMT_bulk,
		Bi214_wire_bulk,
		Bi214_wire_surface,
		Tl208_foil_bulk,
		Tl208_PMT_bulk,
		Pa234m_foil_bulk,
		K40_foil_bulk,
		K40_PMT_bulk,
	]
	for (i, d) in zip(bkgIsotopes, bkgDfs)
		(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( BkgActivityParams, SNparams, SimulationParams, i, d )
		push!(isotope_df, (isotope = string(i), n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	end

	(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( SigActivityParams, SNparams, SimulationParams, :Xi037_foil_bulk, Xi037_foil_bulk )
	push!(isotope_df, (isotope = string(:Xi037_foil_bulk), n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	
	isotope_df
	
	# header = (
	# ["isotope", "expected counts", L"\varepsilon", "activity", "amount", "simulated events"],
	# ["", " ", "[%]", "[Bq/amount]", L"[\textrm{kg or m^3 or l}]", " "]
	# )
	# pretty_table(
	# 	isotope_df,
	# 	header =header,
	# 	backend = Val(:markdown)
	# )
end

# ╔═╡ 94e2f24b-8fb6-4ce3-adfc-532f876a0a5d
#scale histogram to estimated counts
function estimated_counts_hist1D(dataDf, activity, amount, time, binning, nTotalSim)
	totEff = nrow(dataDf) / nTotalSim
	h1d = 
		Hist1D(
			dataDf.reconstructedEnergy1 .+ dataDf.reconstructedEnergy2 ,
			binning
		)
	
	h1d = normalize(h1d) * step(binning) # normalize to area of 1
	
	if(typeof(activity)<:Measurement) # just a precaution, sometime I use activities with uncertainty
		activity = activity.val
	end
	
	h1d.hist.weights = h1d.hist.weights .*(
		activity * time * amount * totEff	
	) # multiply each bin by activity * time * amount * totalEfficiency

	return h1d
end
	

# ╔═╡ cdf251a6-dc5e-4b21-8789-1daafbb676a5
begin 
	# First we normalize histograms to probability
	h1Xi037_foil_bulk = estimated_counts_hist1D( 
		Xi037_foil_bulk, 
		SigActivityParams[:Xi037_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Xi037_foil_bulk] 
	)

	h1Bi214_foil_bulk = estimated_counts_hist1D( 
		Bi214_foil_bulk, 
		BkgActivityParams[:Bi214_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_foil_bulk] 
	)

	h1Bi214_foil_surface = estimated_counts_hist1D( 
		Bi214_foil_surface, 
		BkgActivityParams[:Bi214_foil_surface], 
		SNparams["gasVolume"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_foil_surface] 
	)

	h1Bi214_PMT_bulk = estimated_counts_hist1D( 
		Bi214_PMT_bulk, 
		BkgActivityParams[:Bi214_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_PMT_bulk] 
	)

	h1Bi214_wire_bulk = estimated_counts_hist1D( 
		Bi214_wire_bulk, 
		BkgActivityParams[:Bi214_wire_bulk], # dummy activity
		SNparams["wireBulkMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_wire_bulk] 
	)

	h1Bi214_wire_surface = estimated_counts_hist1D( 
		Bi214_wire_surface, 
		BkgActivityParams[:Bi214_foil_surface], 
		SNparams["gasVolume"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_wire_surface] 
	)

	h1Tl208_foil_bulk = estimated_counts_hist1D( 
		Tl208_foil_bulk, 
		BkgActivityParams[:Tl208_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Tl208_foil_bulk] 
	)

	h1Tl208_PMT_bulk = estimated_counts_hist1D( 
		Tl208_PMT_bulk, 
		BkgActivityParams[:Tl208_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Tl208_PMT_bulk] 
	)

	h1Pa234m_foil_bulk = estimated_counts_hist1D( 
		Pa234m_foil_bulk, 
		BkgActivityParams[:Pa234m_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Pa234m_foil_bulk] 
	)

	h1K40_foil_bulk = estimated_counts_hist1D( 
		K40_foil_bulk, 
		BkgActivityParams[:K40_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:K40_foil_bulk] 
	)

	h1K40_PMT_bulk = estimated_counts_hist1D( 
		K40_PMT_bulk, 
		BkgActivityParams[:K40_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:K40_PMT_bulk] 
	)

	histos = [
		h1Xi037_foil_bulk,
		h1Bi214_foil_bulk,
		h1Bi214_foil_surface,
		h1Bi214_PMT_bulk,
		h1Bi214_wire_bulk,
		h1Bi214_wire_surface,
		h1Tl208_foil_bulk,
		h1Tl208_PMT_bulk,
		h1Pa234m_foil_bulk,
		h1K40_foil_bulk,
		h1K40_PMT_bulk,
	]
	names = [
		"Xi037_foil_bulk",
		"Bi214_foil_bulk",
		"Bi214_foil_surface",
		"Bi214_PMT_bulk",
		"Bi214_wire_bulk",
		"Bi214_wire_surface",
		"Tl208_foil_bulk",
		"Tl208_PMT_bulk",
		"Pa234m_foil_bulk",
		"K40_foil_bulk",
		"K40_PMT_bulk",
	]

	shTotal = SensitivityModule.stackedhist( 
		histos, 
		label = reshape(names, (1, length(names))),
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		title = "Stacked Hist \nEstimated counts in $(SNparams["tYear"]) years of measurement \nusing SDBDRC and recommended activities"
	)

	safesave(plotsdir("SumE", "Stacked_hist_total_spectrum.pdf"), shTotal)
	safesave(plotsdir("SumE", "Stacked_hist_total_spectrum.png"), shTotal)
	
	shTotal
end

# ╔═╡ 8eb7a4e3-5072-43e2-b3c7-d7bdf42eb04b
begin
	p = plot(
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		title = "Estimated counts in $(SNparams["tYear"]) years of measurement \nusing SDBDRC and recommended activities",
		xlims = (0, 4000)
	)
	for i in 1:length(histos)
	plot!(p, 
		histos[i], fa = 0.3, lw = 0, label = names[i]
		#midpoints(binedges(h1Xi037_foil_bulk)), bincounts(h1Xi037_foil_bulk), st = :step, l = names[1] 
	)
	end
	

	safesave(plotsdir("SumE", "hist_total_spectrum.pdf"), p)
	safesave(plotsdir("SumE", "hist_total_spectrum.png"), p)
	
	p
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
	Bi214_foil_bulk_SumE    = 
		Process(Bi214_foil_bulk.reconstructedEnergy1 .+ Bi214_foil_bulk.reconstructedEnergy2, sumEBi214_foil_bulk_Params)
	
	Bi214_foil_surface_SumE = 
		Process(Bi214_foil_surface.reconstructedEnergy1 .+ Bi214_foil_surface.reconstructedEnergy2, sumEBi214_foil_surface_Params)
	
	Bi214_PMT_bulk_SumE     = 
		Process(Bi214_PMT_bulk.reconstructedEnergy1 .+ Bi214_PMT_bulk.reconstructedEnergy2, sumEBi214_PMT_bulk_Params)
	
	Bi214_wire_bulk_SumE    = 
		Process(Bi214_wire_bulk.reconstructedEnergy1 .+ Bi214_wire_bulk.reconstructedEnergy2, sumEBi214_wire_bulk_Params)
	
	Bi214_wire_surface_SumE = 
		Process(Bi214_wire_surface.reconstructedEnergy1 .+ Bi214_wire_surface.reconstructedEnergy2, sumEBi214_wire_surface_Params)
	
	
	Tl208_foil_bulk_SumE    = 
		Process(Tl208_foil_bulk.reconstructedEnergy1 .+ Tl208_foil_bulk.reconstructedEnergy2, sumETl208_foil_bulk_Params)
	
	Tl208_PMT_bulk_SumE     = 
		Process(Tl208_PMT_bulk.reconstructedEnergy1 .+ Tl208_PMT_bulk.reconstructedEnergy2, sumETl208_PMT_bulk_Params)
	
	
	Pa234m_foil_bulk_SumE   = 
		Process(Pa234m_foil_bulk.reconstructedEnergy1 .+ Pa234m_foil_bulk.reconstructedEnergy2, sumEPa234m_foil_bulk_Params)
	
	
	K40_foil_bulk_SumE      = 
		Process(K40_foil_bulk.reconstructedEnergy1 .+ K40_foil_bulk.reconstructedEnergy2, sumEK40_foil_bulk_Params)
	
	K40_PMT_bulk_SumE       = 
		Process(K40_PMT_bulk.reconstructedEnergy1 .+ K40_PMT_bulk.reconstructedEnergy2, sumEK40_PMT_bulk_Params)
	
	
	Xi037_foil_bulk_SumE     = 
		Process(Xi037_foil_bulk.reconstructedEnergy1 .+ Xi037_foil_bulk.reconstructedEnergy2, sumEXiParams)
	

	nothing # so there's no cell output
end

# ╔═╡ ab82a10d-15ca-4528-95e2-e5edfa0dd1c5
md"""
### At construction, the efficiency maps are created as well:
"""

# ╔═╡ d9df2537-7a8a-42dd-a689-6f130d97b525
with(
    gr;

    c = :coolwarm, 
    size = (1600, 600), 
    thickness_scaling = 1.4, 
    right_margin = 8Plots.mm,
    left_margin = 8Plots.mm,
    bottom_margin = 8Plots.mm,
)   do
#Bi214
    p1 = plot(Bi214_foil_bulk_SumE.efficiency, title ="$(Bi214_foil_bulk_SumE.isotopeName)")

    p2 = plot(Bi214_foil_surface_SumE.efficiency, title ="$(Bi214_foil_surface_SumE.isotopeName)")

    p3 = plot(Bi214_wire_bulk_SumE.efficiency, title ="$(Bi214_wire_bulk_SumE.isotopeName)")

    p4 = plot(Bi214_wire_surface_SumE.efficiency, title ="$(Bi214_wire_surface_SumE.isotopeName)")

    p5 = plot(Bi214_PMT_bulk_SumE.efficiency, title ="$(Bi214_PMT_bulk_SumE.isotopeName)")


#Tl208   
    p6 = plot(Tl208_foil_bulk_SumE.efficiency, title ="$(Tl208_foil_bulk_SumE.isotopeName)")

    p7 = plot(Tl208_PMT_bulk_SumE.efficiency, title ="$(Tl208_PMT_bulk_SumE.isotopeName)")

#Pa234m
    p8 = plot(Pa234m_foil_bulk_SumE.efficiency, title ="$(Pa234m_foil_bulk_SumE.isotopeName)")

#K40
    p9 = plot(K40_foil_bulk_SumE.efficiency, title ="$(K40_foil_bulk_SumE.isotopeName)")

    p10 = plot(
		K40_PMT_bulk_SumE.efficiency, 
		title ="$(K40_PMT_bulk_SumE.isotopeName)",
		xlabel = "min sum energy [keV]", 
    	ylabel = "max sum energy [keV]", 
	)

#Xi037
    p11 = plot(
		Xi037_foil_bulk_SumE.efficiency, 
		title ="$(Xi037_foil_bulk_SumE.isotopeName)", 
		xlabel = "min sum energy [keV]", 
    	ylabel = "max sum energy [keV]", 
	)

	plot(p10, p11) # it's not viable to display 11 plots, just look at the saved pngs :P
end


# ╔═╡ ed3a655c-eb27-4cc2-bfb2-6541cd0ca4b5
md"""
#### Combining the **estimated** counts into the efficiency maps we can create the s/b maps:
"""

# ╔═╡ ee5d0838-505b-42e3-89e2-5526de144235
stbSum = get_sToBRatio(
	Xi037_foil_bulk_SumE,
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
);

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


# ╔═╡ 82ec2f54-2044-498f-81e0-9c8c9816d16c
md"""
### The estimated background count in the ROI is: b = $(expBkgESum |> round)
"""

# ╔═╡ df9d6770-a7f0-4e7f-8599-e50703b114c1
effbb = lookup( Xi037_foil_bulk_SumE.efficiency, best_stbSum[:minBinEdge], best_stbSum[:maxBinEdge]);

# ╔═╡ 315e3d2e-a8de-4f56-ac15-0931c3dbfc6f
md"""
### The corresponding signal efficiency in ROI is: $effbb
"""

# ╔═╡ da2d004f-8400-4b90-b4d9-181a1fef318a
md"""
### Finally the calculated sensitivity (assuming Gaussian approximation) for CL = 90% is: 
 ``$T_{1/2}^{2\nu} \geq (3.288 \pm 0.003)\cdot 10^{21} \textrm{yr}$``

"""

# ╔═╡ c5ae4098-ef84-4173-b1d0-d64815e82a0a
ThalfbbESum = get_tHalf(SNparams, effbb, expBkgESum, 1.8)


# ╔═╡ 393c17ab-eb64-4a22-b92b-1b0847a70c1a
begin 
	h1Bkg = sum([
		h1Bi214_foil_bulk,
		h1Bi214_foil_surface,
		h1Bi214_PMT_bulk,
		h1Bi214_wire_bulk,
		h1Bi214_wire_surface,
		h1Tl208_foil_bulk,
		h1Tl208_PMT_bulk,
		h1Pa234m_foil_bulk,
		h1K40_foil_bulk,
		h1K40_PMT_bulk,
	])

	SensitivityModule.stackedstephist(
		[h1Xi037_foil_bulk, h1Bkg], 
		label=["signal" "background"], 
		xlabel= L"E_{sum} "*"[keV]", 
		ylabel = "counts [$(step(binedges(h1Xi037_foil_bulk))) "* L"\textrm{keV ^{-1}}]", 
		title="stacked histogram signal and background",
		c=[1 3]
	)

	safesave(plotsdir("SumE","h1_step_bkg_sig.png"), current())
	current()

	h1BkgROI = restrict(h1Bkg, 800, 1500)
	h1SigROI = restrict(h1Xi037_foil_bulk, 800, 1500)
	
	bkgCtsROI = round(integral( h1BkgROI ), digits = 2)
	sigEffROI = round(lookup( Xi037_foil_bulk_SumE.efficiency, 800, 1500 )*100, sigdigits=3) 
	
	SensitivityModule.stackedhist!([h1SigROI, h1BkgROI], c=[1 3], label= ["sig ∈ (800, 1500)keV" "bkg ∈ (800, 1500)keV"])

	annotate!([
		(
			2050, 1.25e4, "signal ε = $sigEffROI%", :left
		),
		(
			2050, 1.15e4, "expected b = $bkgCtsROI", :left
		),
	])

	safesave(plotsdir("SumE","h1_step_bkg_sig_ROI.png"), current())
	current()
	
end

# ╔═╡ be2717cc-0a28-4012-86df-f419a46f629b


# ╔═╡ 70b727b5-a0f0-4d4b-903a-1d61eb433ee4
length(midpoints(binedges(h1Xi037_foil_bulk)))

# ╔═╡ Cell order:
# ╠═af92ba46-999c-11ee-11a4-8554ebc58c0d
# ╠═fcfd9507-542c-4bab-9b8b-5e8a00b6bb55
# ╠═8678f310-e581-41ec-92f2-d5b4e1461396
# ╠═cf56ae15-6048-48fb-ae49-a8c7cd7390eb
# ╟─8d36f771-6342-4237-a3fc-b10727fc18e0
# ╠═f93ae5c5-0702-4ec2-9294-7d506153c81e
# ╠═3a279499-3b4c-4636-b9a6-6183868c3ab9
# ╠═d9102629-c1d6-4a05-aae4-af49d9c9ddbb
# ╠═ed904b40-d638-44a3-8e01-1bd86dd1a919
# ╟─86932308-47bf-4ec9-bd53-8cc6b132d32b
# ╠═e1015a40-3307-46b4-8099-1cac469350a3
# ╟─86db346d-5e02-444e-aa1a-892e64bc42b3
# ╟─01fdc8b3-e467-401e-8ccb-8397b72d057a
# ╟─a5f981b8-aa70-4d1d-a82b-491feb1a1017
# ╠═89852571-c350-4f50-84e0-5eeba29b5da7
# ╠═94e2f24b-8fb6-4ce3-adfc-532f876a0a5d
# ╠═cdf251a6-dc5e-4b21-8789-1daafbb676a5
# ╠═8eb7a4e3-5072-43e2-b3c7-d7bdf42eb04b
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
# ╠═da2d004f-8400-4b90-b4d9-181a1fef318a
# ╠═c5ae4098-ef84-4173-b1d0-d64815e82a0a
# ╠═393c17ab-eb64-4a22-b92b-1b0847a70c1a
# ╠═be2717cc-0a28-4012-86df-f419a46f629b
# ╠═70b727b5-a0f0-4d4b-903a-1d61eb433ee4
