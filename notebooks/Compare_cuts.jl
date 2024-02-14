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
	include(scriptsdir("LoadData.jl"))
end

# ╔═╡ 71739742-f9e8-47c1-8f16-326f73bce8c2
using ColorSchemes

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

# ╔═╡ 0f712dc8-2725-447d-be51-82135b04ac9a
md"""
# Isotope details using SDBDRC
"""

# ╔═╡ e4ed8dd5-645d-4e39-8efe-a9b9d1f17376
md"""
# Isotope details using SDBDRC + vertex
"""

# ╔═╡ 23213a05-cb68-4f16-be98-510cfede5bcd
begin
	isotope_vertex_df = DataFrame()
	for key in keys(vertexFiles)
		if ( occursin("bb", key) || occursin("Xi", key))
			(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( 
				SigActivityParams, SNparams, SimulationParams, Symbol(key), vertexFiles[key] 
			)
		else
			(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( 
				BkgActivityParams, SNparams, SimulationParams, Symbol(key), vertexFiles[key] 
			)
		end
		push!(isotope_vertex_df, (isotope = key, n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	end

	push!(
		isotope_vertex_df, 
		(
			isotope = "Total_bkg", 
			n_expected = sum(isotope_vertex_df[1:end-1,2]), 
			ε = 0.0, 
			activity = 0.0, 
			amount = 0.0, 
			n_simulated = 0.0
		), 
		promote=:true
	)
	
	# isotope_vertex_df
	

end

# ╔═╡ b0895c4b-bda5-46bd-b3bf-ecfe4b806353
begin
	isotope_SDBDRC_df = DataFrame()
	for key in keys(SDBDRCFiles)
		if ( occursin("bb", key) || occursin("Xi", key))
			(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( 
				SigActivityParams, SNparams, SimulationParams, Symbol(key), SDBDRCFiles[key] 
			)
		else
			(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( 
				BkgActivityParams, SNparams, SimulationParams, Symbol(key), SDBDRCFiles[key] 
			)
		end
		push!(isotope_SDBDRC_df, (isotope = key, n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	end

	push!(
		isotope_SDBDRC_df, 
		(
			isotope = "Total_bkg", 
			n_expected = sum(isotope_SDBDRC_df[1:end-1,2]), 
			ε = 0.0, 
			activity = 0.0, 
			amount = 0.0, 
			n_simulated = 0.0
		), 
		promote=:true
	)
	
	isotope_SDBDRC_df
	
	header = (
	["isotope", "expected counts", L"\varepsilon", "activity", "amount", "simulated events"],
	["", " ", "[%]", "[Bq/amount]", L"[\textrm{kg or m^3 or l}]", " "]
	)
	pretty_table(
		isotope_vertex_df,
		header =header,
		backend = Val(:markdown)
	)
end

# ╔═╡ 24266b27-e3c2-4640-9a9e-447045364d81
md"""
# Isotope details using SDBDRC + vertex + probability
"""

# ╔═╡ 1225796b-2d3a-4842-8c1a-89718811db52
begin
	isotope_prob_df = DataFrame()
	for key in keys(probFiles)
		if ( occursin("bb", key) || occursin("Xi", key))
			(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( 
				SigActivityParams, SNparams, SimulationParams, Symbol(key), probFiles[key] 
			)
		else
			(a, m, t, ε, nExpTot, nTotSim) = get_isotope_details( 
				BkgActivityParams, SNparams, SimulationParams, Symbol(key), probFiles[key] 
			)
		end
		push!(isotope_prob_df, (isotope = key, n_expected = Measurements.value(nExpTot), ε = ε*100, activity = Measurements.value(a), amount = m, n_simulated = nTotSim), promote=:true)
	end

	push!(
		isotope_prob_df, 
		(
			isotope = "Total_bkg", 
			n_expected = sum(isotope_prob_df[1:end-1,2]), 
			ε = 0.0, 
			activity = 0.0, 
			amount = 0.0, 
			n_simulated = 0.0
		), 
		promote=:true
	)
	
	isotope_prob_df

	# pretty_table(
	# 	isotope_prob_df,
	# 	header =header,
	# 	backend = Val(:markdown)
	# )
end

# ╔═╡ 4a1b3fdf-1a88-4fc4-b670-e62db0d85573
md"""
# Estimated counts Hist1D
"""

# ╔═╡ 6fcbfdb0-3496-4edd-a378-cb02828c84d5
SDBDRCFiles["Xi037_foil_bulk"]

# ╔═╡ 078b195e-e83d-4898-a593-bc704fe3bb9b
function get_estimated_h1( data, nExp, binning )
	h = Hist1D(data, binning)
	h = normalize(h, width =:false)

	h.hist.weights .*= nExp
	return h
end

# ╔═╡ 5efd324d-8fdf-4d17-a09a-cbd4ea277312
cScheme = ColorSchemes.tab20

# ╔═╡ 55be8367-6822-4ad2-9c6f-cd66327dc1c9
begin 
	histosSDBDRC = []
	labelsSDBDRC = []
	h1sigSDBDRC = nothing
	sigLabel = nothing
	
	for (r, key) in enumerate(keys(SDBDRCFiles))
		data = SDBDRCFiles[key].reconstructedEnergy1 .+ SDBDRCFiles[key].reconstructedEnergy2
		if ( occursin("bb", key))
			continue
		end
		if (occursin("Xi", key))
			h1sigSDBDRC = get_estimated_h1( data, isotope_SDBDRC_df[r ,2], sumEParams[:binning]  )
			sigLabel = key
			continue
		end
		hh = get_estimated_h1( data, isotope_SDBDRC_df[r ,2], sumEParams[:binning]  )
		push!(histosSDBDRC, hh)
		push!(labelsSDBDRC, key)
	end

	labelsSDBDRC = vcat(sigLabel, labelsSDBDRC)
	shTotalSDBDRC = SensitivityModule.stackedhist( 
		vcat(h1sigSDBDRC, histosSDBDRC ), 
		label = reshape(labelsSDBDRC, (1, length(labelsSDBDRC))),
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		title = "stacked hist estimated counts \n $(SNparams["tYear"]) yr with recommended activities \nusing SDBDRC cuts",
		c = [cScheme[i] for i=1:length(vcat(h1sigSDBDRC, histosSDBDRC ))]',
		yformatter = :scientific,
		ylims = (0, 2.2e4),
		legend = :outerright,
		size = (2200, 1000),
		thickness_scaling = 1.8
	)
	safesave( plotsdir("CompareCuts", "SDBDRCcuts.pdf"), current() )
	safesave( plotsdir("CompareCuts", "SDBDRCcuts.svg"), current() )
	current()
end


# ╔═╡ 467bc8f3-f65d-42c8-9ce2-156adbd8f670
begin 
	histosvertex = []
	labelsvertex = []
	h1sigvertex = nothing
	
	for (r, key) in enumerate(keys(vertexFiles))
		data = vertexFiles[key].reconstructedEnergy1 .+ vertexFiles[key].reconstructedEnergy2
		if ( occursin("bb", key))
			continue
		end
		if (occursin("Xi", key))
			h1sigvertex = get_estimated_h1( data, isotope_vertex_df[r ,2], sumEParams[:binning]  )
			continue
		end
		hh = get_estimated_h1( data, isotope_vertex_df[r ,2], sumEParams[:binning]  )
		push!(histosvertex, hh)
		push!(labelsvertex, key)
	end

	labelsvertex = vcat(sigLabel, labelsvertex)
	shTotalvertex = SensitivityModule.stackedhist( 
		vcat(h1sigvertex, histosvertex), 
		label = reshape(labelsvertex, (1, length(labelsvertex))),
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		title = "stacked hist estimated counts \n $(SNparams["tYear"]) yr with recommended activities \nusing SDBDRC + 100mm vertex cuts",
		c = [cScheme[i] for i=1:length(labelsvertex)]',
		yformatter = :scientific,
		ylims = (0, 2.2e4),
		legend = :outerright,
		size = (2200, 1000),
		thickness_scaling = 1.8
	)

	safesave( plotsdir("CompareCuts", "vertexcuts.pdf"), current() )
	safesave( plotsdir("CompareCuts", "vertexcuts.svg"), current() )
	current()
end

# ╔═╡ d477a745-e791-4706-971a-19d22958c158
begin 
	histosprob = []
	labelsprob = []
	h1sigprob = nothing

	
	for (r, key) in enumerate(keys(probFiles))
		data = probFiles[key].reconstructedEnergy1 .+ probFiles[key].reconstructedEnergy2
		if ( occursin("bb", key))
			continue
		end
		if (occursin("Xi", key))
			h1sigprob = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
			continue
		end
		hh = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
		push!(histosprob, hh)
		push!(labelsprob, key)
	end
	
	labelsbkg = vcat(sigLabel, labelsprob)
	
	shTotalprob = SensitivityModule.stackedhist( 
		vcat(h1sigprob, histosprob), 
		label = reshape(labelsbkg, (1, length(labelsbkg))),
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		yformatter = :scientific,
		title = "stacked hist estimated counts \n $(SNparams["tYear"]) yr with recommended activities \nusing SDBDRC + 100mm vertex + probability cuts",
		c = [cScheme[i] for i=1:length(labelsbkg)]',
		ylims = (0, 2.2e4),
		# ylims = (0, 1e4),
		legend = :outerright,
		size = (2200, 1000),
		thickness_scaling = 1.8
	)

	safesave( plotsdir("CompareCuts", "probcuts.pdf"), current() )
	safesave( plotsdir("CompareCuts", "probcuts.svg"), current() )
	current()
end

# ╔═╡ e28ebb26-29b4-45c1-a0db-92d63e01500e
begin 
	histosprobbg = []
	labelsprobbg = []
	histosprobsg = []
	labelsprobsg = []
	for (r, key) in enumerate(keys(probFiles))
		data = probFiles[key].reconstructedEnergy1 .+ probFiles[key].reconstructedEnergy2
		if ( occursin("bb", key))
			continue
		end
		if ( occursin("Xi", key))
			hh = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
			push!(histosprobsg, hh)
			push!(labelsprobsg, key)
			continue
		end
		hh = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
		push!(histosprobbg, hh)
		push!(labelsprobsg, key)
	end
	
	hhbg= sum(histosprobbg)
	
	@show integral(histosprobsg[1], width = true)
	
	shTotalprob1 = SensitivityModule.stackedhist( 
		[histosprobsg[1], hhbg], 
		label = ["Xi31" "sum of backgrounds"],
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		title = "Stacked Hist \nEstimated counts in $(SNparams["tYear"]) years of measurement \nusing SDBDRC + vertex + prob and \nrecommended activities"
	)
end

# ╔═╡ 8adb3f23-19e2-440d-a90d-3ca81a7ea712
begin 
	PMT = Hist1D(Float64; bins = sumEParams[:binning])
	foil_bulk = Hist1D(Float64; bins = sumEParams[:binning])
	radon = Hist1D(Float64; bins = sumEParams[:binning])
	labelsPMT = []
	labelsfoil = []
	labelsradon = []
	
	for (r, key) in enumerate(keys(probFiles))
		data = probFiles[key].reconstructedEnergy1 .+ probFiles[key].reconstructedEnergy2
		if ( occursin("bb", key))
			continue
		end
		if ( occursin("Xi", key))
			continue
		end
		if (occursin("PMT", key))
			hh = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
			PMT += hh
			push!(labelsPMT, key)
			continue
		end
		if (occursin("foil_bulk", key))
			hh = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
			foil_bulk += hh
			push!(labelsfoil, key)
			continue
		end
		hh = get_estimated_h1( data, isotope_prob_df[r ,2], sumEParams[:binning]  )
		radon += hh
		push!(labelsradon, key)
	end

	nS =  integral(h1sigprob) |> round 
	nR =  integral(radon)  |> round
	nP =  integral(PMT) |> round
	nF =  integral(foil_bulk) |> round
	
	labs = ["signal; $nS" "PMT bkg; $nP" "foil_bulk bkg; $nF" "radon bkg; $nR"]



	SensitivityModule.stackedhist(
		vcat(h1sigprob, PMT, foil_bulk, radon),
		label = labs,
		legendtitle ="source; entries",
		legendtitlefontsize = 16,
		legendfontsize = 16,
		xlabel = L"E_{sum}"*" [keV]",
		ylabel = "estimated counts " *L"[%$(step(sumEParams[:binning])) \textrm{keV^{-1}}]",
		yformatter = :scientific,
		title = "stacked hist estimated counts \n $(SNparams["tYear"]) yr with recommended activities \nusing SDBDRC + 100mm vertex + probability cuts",
		c = [1, 2, 10, 11]',
		ylims = (0, 1e4),
		legend = :outerright,
		size = (2200, 1000),
		thickness_scaling = 1.8
	)

	safesave( plotsdir("CompareCuts", "probcuts_srcs.pdf"), current() )
	safesave( plotsdir("CompareCuts", "probcuts_srcs.svg"), current() )
	current()
end


# ╔═╡ f80a0885-1ce0-4736-b1ac-e2fd7854a8a5
savefig

# ╔═╡ Cell order:
# ╠═60342f92-c102-11ee-04ee-4fe6bdbcb411
# ╠═a710d50c-7c23-4621-8a25-41ba31896c16
# ╠═a902895b-0b5e-45c7-8eb2-c05c472556cb
# ╠═79a8065f-fc6b-49dd-a3be-b21fb343b6f7
# ╠═c9c72836-2825-4e9b-ba28-57d102ec8b33
# ╟─0f712dc8-2725-447d-be51-82135b04ac9a
# ╠═b0895c4b-bda5-46bd-b3bf-ecfe4b806353
# ╠═e4ed8dd5-645d-4e39-8efe-a9b9d1f17376
# ╠═23213a05-cb68-4f16-be98-510cfede5bcd
# ╠═24266b27-e3c2-4640-9a9e-447045364d81
# ╠═1225796b-2d3a-4842-8c1a-89718811db52
# ╟─4a1b3fdf-1a88-4fc4-b670-e62db0d85573
# ╠═6fcbfdb0-3496-4edd-a378-cb02828c84d5
# ╠═078b195e-e83d-4898-a593-bc704fe3bb9b
# ╠═71739742-f9e8-47c1-8f16-326f73bce8c2
# ╠═5efd324d-8fdf-4d17-a09a-cbd4ea277312
# ╠═55be8367-6822-4ad2-9c6f-cd66327dc1c9
# ╠═467bc8f3-f65d-42c8-9ce2-156adbd8f670
# ╠═d477a745-e791-4706-971a-19d22958c158
# ╠═e28ebb26-29b4-45c1-a0db-92d63e01500e
# ╠═8adb3f23-19e2-440d-a90d-3ca81a7ea712
# ╠═f80a0885-1ce0-4736-b1ac-e2fd7854a8a5
