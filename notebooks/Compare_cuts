### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 60342f92-c102-11ee-04ee-4fe6bdbcb411
using DrWatson

# ╔═╡ a710d50c-7c23-4621-8a25-41ba31896c16
@quickactivate "SNSensitivityEstimate"

# ╔═╡ a902895b-0b5e-45c7-8eb2-c05c472556cb
begin
	using StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Measurements, PrettyTables
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
end

# ╔═╡ 79a8065f-fc6b-49dd-a3be-b21fb343b6f7
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

# ╔═╡ c9c72836-2825-4e9b-ba28-57d102ec8b33
ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"] ) 

# ╔═╡ 8a70b916-12ba-43a9-b834-fb07fa7b4339
begin
#Bi214
	Bi214_foil_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Bi214_foil_bulk_1e8E.root")) |> ffrf
	Bi214_foil_surface_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Bi214_foil_surface_1e7E.root")) |> ffrf
	Bi214_PMT_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Bi214_PMT_bulk_1e7E.root")) |> ffrf
	Bi214_wire_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Bi214_wire_bulk_1e7E.root")) |> ffrf
	Bi214_wire_surface_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Bi214_wire_surface_1e7E.root")) |> ffrf

#Tl208
	Tl208_foil_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Tl208_foil_bulk_1e8E.root")) |> ffrf
	Tl208_PMT_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Tl208_PMT_bulk_1e8E.root")) |> ffrf

#Pa234m
	Pa234m_foil_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Pa234m_foil_bulk_1e8E.root")) |> ffrf

#K40
	K40_foil_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/K40_foil_bulk_1e8E.root")) |> ffrf
	K40_PMT_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/K40_PMT_bulk_1e8E.root")) |> ffrf

#Xi31
	Xi037_foil_bulk_SDBDRC = 
		ROOTFile(datadir("sims/SDBDRC/Xi037_foil_bulk_1e8E.root")) |> ffrf 
end

# ╔═╡ fef51ff8-f608-4976-a7d6-0215c94fb0bd
begin
#Bi214
	Bi214_foil_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Bi214_foil_bulk_1e8E_SDBDRC_vertex.root")) |> ffrf
	Bi214_foil_surface_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Bi214_foil_surface_1e7E_SDBDRC_vertex.root")) |> ffrf
	Bi214_PMT_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Bi214_PMT_bulk_1e7E_SDBDRC_vertex.root")) |> ffrf
	Bi214_wire_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Bi214_wire_bulk_1e7E_SDBDRC_vertex.root")) |> ffrf
	Bi214_wire_surface_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Bi214_wire_surface_1e7E_SDBDRC_vertex.root")) |> ffrf

#Tl208
	Tl208_foil_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Tl208_foil_bulk_1e8E_SDBDRC_vertex.root")) |> ffrf
	Tl208_PMT_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Tl208_PMT_bulk_1e8E_SDBDRC_vertex.root")) |> ffrf

#Pa234m
	Pa234m_foil_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/Pa234m_foil_bulk_1e8E_SDBDRC_vertex.root")) |> ffrf

#K40
	K40_foil_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/K40_foil_bulk_1e8E_SDBDRC_vertex.root")) |> ffrf
	K40_PMT_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC_vertex/K40_PMT_bulk_1e8E_SDBDRC_vertex.root")) |> ffrf

#Xi31
	Xi037_foil_bulk_SDBDRC_vertex = 
		ROOTFile(datadir("sims/SDBDRC/Xi037_foil_bulk_1e8E.root")) |> ffrf # NOT YET SIMULATED, USING OLD FILE!!!
end

# ╔═╡ 0f712dc8-2725-447d-be51-82135b04ac9a
md"""
# Isotope details using SDBDRC
"""

# ╔═╡ b0895c4b-bda5-46bd-b3bf-ecfe4b806353
let
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
		Bi214_foil_bulk_SDBDRC,
		Bi214_foil_surface_SDBDRC,
		Bi214_PMT_bulk_SDBDRC,
		Bi214_wire_bulk_SDBDRC,
		Bi214_wire_surface_SDBDRC,
		Tl208_foil_bulk_SDBDRC,
		Tl208_PMT_bulk_SDBDRC,
		Pa234m_foil_bulk_SDBDRC,
		K40_foil_bulk_SDBDRC,
		K40_PMT_bulk_SDBDRC,
	]
	for (i, d) in zip(bkgIsotopes, bkgDfs)
		(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( BkgActivityParams, SNparams, SimulationParams, i, d )
		push!(isotope_df, (isotope = string(i), n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	end

	(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( SigActivityParams, SNparams, SimulationParams, :Xi037_foil_bulk, Xi037_foil_bulk_SDBDRC )
	
	push!(isotope_df, (isotope = string(:Xi037_foil_bulk), n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)

	push!(isotope_df, (
	isotope = "Total_bkg", 
	n_expected = sum(isotope_df[1:end-1,2]), 
	ε = 0.0, 
	activity = 0.0, 
	amount = 0.0, 
	n_simulated = 0.0
	), promote=:true)
	
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

# ╔═╡ e4ed8dd5-645d-4e39-8efe-a9b9d1f17376
md"""
# Isotope details using SDBDRC + vertex
"""

# ╔═╡ 23213a05-cb68-4f16-be98-510cfede5bcd
let
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
		Bi214_foil_bulk_SDBDRC_vertex,
		Bi214_foil_surface_SDBDRC_vertex,
		Bi214_PMT_bulk_SDBDRC_vertex,
		Bi214_wire_bulk_SDBDRC_vertex,
		Bi214_wire_surface_SDBDRC_vertex,
		Tl208_foil_bulk_SDBDRC_vertex,
		Tl208_PMT_bulk_SDBDRC_vertex,
		Pa234m_foil_bulk_SDBDRC_vertex,
		K40_foil_bulk_SDBDRC_vertex,
		K40_PMT_bulk_SDBDRC_vertex,
	]
	for (i, d) in zip(bkgIsotopes, bkgDfs)
		(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( BkgActivityParams, SNparams, SimulationParams, i, d )
		push!(isotope_df, (isotope = string(i), n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	end

	(a, m, t, ε, nExpTot, nTotSim) = 
		get_isotope_details( 
			SigActivityParams, SNparams, SimulationParams, :Xi037_foil_bulk, Xi037_foil_bulk_SDBDRC_vertex 
		)
	
	push!(
		isotope_df, 
		(
			isotope = string(:Xi037_foil_bulk), 
			n_expected = Measurements.value(nExpTot), 
			ε = ε*100, 
			activity = Measurements.value(a), 
			amount = m, 
			n_simulated = nTotSim
		), 
		promote=:true)

	push!(isotope_df, (
		isotope = "Total_bkg", 
		n_expected = sum(isotope_df[1:end-1,2]), 
		ε = 0.0, 
		activity = 0.0, 
		amount = 0.0, 
		n_simulated = 0.0
		), promote=:true
	)
	
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

# ╔═╡ 4a1b3fdf-1a88-4fc4-b670-e62db0d85573
md"""
# Estimated counts Hist1D
"""

# ╔═╡ 467bc8f3-f65d-42c8-9ce2-156adbd8f670
begin 
	# First we normalize histograms to probability
	h1Xi037_foil_bulk_SDBDRC = estimated_counts_hist1D( 
		Xi037_foil_bulk_SDBDRC, 
		SigActivityParams[:Xi037_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Xi037_foil_bulk] 
	)

	h1Bi214_foil_bulk_SDBDRC = estimated_counts_hist1D( 
		Bi214_foil_bulk_SDBDRC, 
		BkgActivityParams[:Bi214_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_foil_bulk] 
	)

	h1Bi214_foil_surface_SDBDRC = estimated_counts_hist1D( 
		Bi214_foil_surface_SDBDRC, 
		BkgActivityParams[:Bi214_foil_surface], 
		SNparams["gasVolume"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_foil_surface] 
	)

	h1Bi214_PMT_bulk_SDBDRC = estimated_counts_hist1D( 
		Bi214_PMT_bulk_SDBDRC, 
		BkgActivityParams[:Bi214_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_PMT_bulk] 
	)

	h1Bi214_wire_bulk_SDBDRC = estimated_counts_hist1D( 
		Bi214_wire_bulk_SDBDRC, 
		BkgActivityParams[:Bi214_wire_bulk], # dummy activity
		SNparams["wireBulkMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_wire_bulk] 
	)

	h1Bi214_wire_surface_SDBDRC = estimated_counts_hist1D( 
		Bi214_wire_surface_SDBDRC, 
		BkgActivityParams[:Bi214_foil_surface], 
		SNparams["gasVolume"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_wire_surface] 
	)

	h1Tl208_foil_bulk_SDBDRC = estimated_counts_hist1D( 
		Tl208_foil_bulk_SDBDRC, 
		BkgActivityParams[:Tl208_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Tl208_foil_bulk] 
	)

	h1Tl208_PMT_bulk_SDBDRC = estimated_counts_hist1D( 
		Tl208_PMT_bulk_SDBDRC, 
		BkgActivityParams[:Tl208_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Tl208_PMT_bulk] 
	)

	h1Pa234m_foil_bulk_SDBDRC = estimated_counts_hist1D( 
		Pa234m_foil_bulk_SDBDRC, 
		BkgActivityParams[:Pa234m_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Pa234m_foil_bulk] 
	)

	h1K40_foil_bulk_SDBDRC = estimated_counts_hist1D( 
		K40_foil_bulk_SDBDRC, 
		BkgActivityParams[:K40_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:K40_foil_bulk] 
	)

	h1K40_PMT_bulk_SDBDRC = estimated_counts_hist1D( 
		K40_PMT_bulk_SDBDRC, 
		BkgActivityParams[:K40_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:K40_PMT_bulk] 
	)

	histos = [
		h1Xi037_foil_bulk_SDBDRC,
		h1Bi214_foil_bulk_SDBDRC,
		h1Bi214_foil_surface_SDBDRC,
		h1Bi214_PMT_bulk_SDBDRC,
		h1Bi214_wire_bulk_SDBDRC,
		h1Bi214_wire_surface_SDBDRC,
		h1Tl208_foil_bulk_SDBDRC,
		h1Tl208_PMT_bulk_SDBDRC,
		h1Pa234m_foil_bulk_SDBDRC,
		h1K40_foil_bulk_SDBDRC,
		h1K40_PMT_bulk_SDBDRC,
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

# ╔═╡ d477a745-e791-4706-971a-19d22958c158
begin 
	# First we normalize histograms to probability
	h1Xi037_foil_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Xi037_foil_bulk_SDBDRC_vertex, 
		SigActivityParams[:Xi037_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Xi037_foil_bulk] 
	)

	h1Bi214_foil_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Bi214_foil_bulk_SDBDRC_vertex, 
		BkgActivityParams[:Bi214_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_foil_bulk] 
	)

	h1Bi214_foil_surface_SDBDRC_vertex = estimated_counts_hist1D( 
		Bi214_foil_surface_SDBDRC_vertex, 
		BkgActivityParams[:Bi214_foil_surface], 
		SNparams["gasVolume"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_foil_surface] 
	)

	h1Bi214_PMT_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Bi214_PMT_bulk_SDBDRC_vertex, 
		BkgActivityParams[:Bi214_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_PMT_bulk] 
	)

	h1Bi214_wire_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Bi214_wire_bulk_SDBDRC_vertex, 
		BkgActivityParams[:Bi214_wire_bulk], # dummy activity
		SNparams["wireBulkMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_wire_bulk] 
	)

	h1Bi214_wire_surface_SDBDRC_vertex = estimated_counts_hist1D( 
		Bi214_wire_surface_SDBDRC_vertex, 
		BkgActivityParams[:Bi214_foil_surface], 
		SNparams["gasVolume"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Bi214_wire_surface] 
	)

	h1Tl208_foil_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Tl208_foil_bulk_SDBDRC_vertex, 
		BkgActivityParams[:Tl208_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Tl208_foil_bulk] 
	)

	h1Tl208_PMT_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Tl208_PMT_bulk_SDBDRC_vertex, 
		BkgActivityParams[:Tl208_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Tl208_PMT_bulk] 
	)

	h1Pa234m_foil_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		Pa234m_foil_bulk_SDBDRC_vertex, 
		BkgActivityParams[:Pa234m_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:Pa234m_foil_bulk] 
	)

	h1K40_foil_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		K40_foil_bulk_SDBDRC_vertex, 
		BkgActivityParams[:K40_foil_bulk], 
		SNparams["foilMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:K40_foil_bulk] 
	)

	h1K40_PMT_bulk_SDBDRC_vertex = estimated_counts_hist1D( 
		K40_PMT_bulk_SDBDRC_vertex, 
		BkgActivityParams[:K40_PMT_bulk], 
		SNparams["PMTGlassMass"], 
		SNparams["t"], 
		sumEParams[:binning], 
		SimulationParams[:K40_PMT_bulk] 
	)

	histos_vertex = [
		h1Xi037_foil_bulk_SDBDRC_vertex,
		h1Bi214_foil_bulk_SDBDRC_vertex,
		h1Bi214_foil_surface_SDBDRC_vertex,
		h1Bi214_PMT_bulk_SDBDRC_vertex,
		h1Bi214_wire_bulk_SDBDRC_vertex,
		h1Bi214_wire_surface_SDBDRC_vertex,
		h1Tl208_foil_bulk_SDBDRC_vertex,
		h1Tl208_PMT_bulk_SDBDRC_vertex,
		h1Pa234m_foil_bulk_SDBDRC_vertex,
		h1K40_foil_bulk_SDBDRC_vertex,
		h1K40_PMT_bulk_SDBDRC_vertex,
	]

	shTotal_vertex = SensitivityModule.stackedhist( 
		histos, 
		label = reshape(names, (1, length(names))),
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		title = "Stacked Hist \nEstimated counts in $(SNparams["tYear"]) years of measurement \nusing SDBDRC + vertex cut and recommended activities"
	)

	safesave(plotsdir("SumE", "Stacked_hist_total_spectrum_vertex.pdf"), shTotal)
	safesave(plotsdir("SumE", "Stacked_hist_total_spectrum_vertex.png"), shTotal)
	
	shTotal_vertex
end

# ╔═╡ 8cf24c42-3333-4004-9b6d-85d0c67903d3
begin 
	h1Bkg_SDBDRC = sum([
		h1Bi214_foil_bulk_SDBDRC,
		h1Bi214_foil_surface_SDBDRC,
		h1Bi214_PMT_bulk_SDBDRC,
		h1Bi214_wire_bulk_SDBDRC,
		h1Bi214_wire_surface_SDBDRC,
		h1Tl208_foil_bulk_SDBDRC,
		h1Tl208_PMT_bulk_SDBDRC,
		h1Pa234m_foil_bulk_SDBDRC,
		h1K40_foil_bulk_SDBDRC,
		h1K40_PMT_bulk_SDBDRC,
	])

	SensitivityModule.stackedstephist(
		[h1Xi037_foil_bulk_SDBDRC, h1Bkg_SDBDRC], 
		label=["signal" "background"], 
		xlabel= L"E_{sum} "*"[keV]", 
		ylabel = "counts [$(step(binedges(h1Xi037_foil_bulk_SDBDRC))) "* L"\textrm{keV ^{-1}}]", 
		title="stacked histogram signal and background \nSDBDRC cuts",
		c=[1 3]
	)

	h1Bkg_SDBDRC_vertex = sum([
		h1Bi214_foil_bulk_SDBDRC_vertex,
		h1Bi214_foil_surface_SDBDRC_vertex,
		h1Bi214_PMT_bulk_SDBDRC_vertex,
		h1Bi214_wire_bulk_SDBDRC_vertex,
		h1Bi214_wire_surface_SDBDRC_vertex,
		h1Tl208_foil_bulk_SDBDRC_vertex,
		h1Tl208_PMT_bulk_SDBDRC_vertex,
		h1Pa234m_foil_bulk_SDBDRC_vertex,
		h1K40_foil_bulk_SDBDRC_vertex,
		h1K40_PMT_bulk_SDBDRC_vertex,
	])

	SensitivityModule.stackedstephist!(
		[h1Xi037_foil_bulk_SDBDRC_vertex, h1Bkg_SDBDRC_vertex], 
		label=["signal + vertex" "background + vertex"], 
		xlabel= L"E_{sum} "*"[keV]", 
		ylabel = "counts [$(step(binedges(h1Xi037_foil_bulk_SDBDRC_vertex))) "* L"\textrm{keV ^{-1}}]", 
		title="stacked histogram signal and background \nSDBDRC cuts + vertex",
		c=[2 4]
	)

	safesave(plotsdir("SumE","h1_step_bkg_sig_SDBDRC_vertex_compare.png"), current())
	current()
	
end

# ╔═╡ Cell order:
# ╠═60342f92-c102-11ee-04ee-4fe6bdbcb411
# ╠═a710d50c-7c23-4621-8a25-41ba31896c16
# ╠═a902895b-0b5e-45c7-8eb2-c05c472556cb
# ╠═79a8065f-fc6b-49dd-a3be-b21fb343b6f7
# ╠═c9c72836-2825-4e9b-ba28-57d102ec8b33
# ╠═8a70b916-12ba-43a9-b834-fb07fa7b4339
# ╠═fef51ff8-f608-4976-a7d6-0215c94fb0bd
# ╟─0f712dc8-2725-447d-be51-82135b04ac9a
# ╟─b0895c4b-bda5-46bd-b3bf-ecfe4b806353
# ╟─e4ed8dd5-645d-4e39-8efe-a9b9d1f17376
# ╟─23213a05-cb68-4f16-be98-510cfede5bcd
# ╟─4a1b3fdf-1a88-4fc4-b670-e62db0d85573
# ╠═467bc8f3-f65d-42c8-9ce2-156adbd8f670
# ╠═d477a745-e791-4706-971a-19d22958c158
# ╠═8cf24c42-3333-4004-9b6d-85d0c67903d3
