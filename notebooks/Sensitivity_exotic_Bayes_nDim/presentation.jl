### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 32c435b2-ee35-4a65-ba59-b23ee5ba7992
using DrWatson

# ╔═╡ 432c7c60-37ba-4ade-8935-500652ddc0c3
begin
	# using DrWatson
	@quickactivate "SNSensitivityEstimate"
	using CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Measurements, PrettyTables, BAT, PlutoUI, LaTeXStrings
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
	nothing
end

# ╔═╡ 5479bdcf-bb42-4910-b061-befdb1b86b37
html"<button onclick='present()'>present</button>"

# ╔═╡ 64580852-4944-45cb-9b9f-b087edb1c7ff
html"""
<style>
  main {
    max-width: 800px;
	font-size: 18px;
	}
</style>
"""

# ╔═╡ 0bd50ff5-d407-4ba1-b4ab-f5fc52ee3c17
md"""
# Updates on the sensitivity study to BSM physics.
### New and improved multi-dimensional frequentist limits & improved Bayesian models

#### Maros Petro
"""

# ╔═╡ 137994dd-199b-4531-98f0-ba29a44bc4db
md"""
# Outline

### 1. Background model
### 2.  Methodolofy
   1. #### Feldman Cousins sensitivity
   2. #### Bayesian sensitivity
### 3. Results
   1. #### $$0\nu\beta\beta$$
   2. #### $$0\nu\beta\beta\chi^0$$ & $$0\nu\beta\beta\chi^0\chi^0$$
   3. #### $$\nu_R\nu_L\beta\beta$$
   4. #### $$2\nu\beta\beta$$ with $$\xi_{31}, \xi_{51}$$
"""

# ╔═╡ cea77a5a-92f3-412d-8e50-69f4252eb52a
md"""
# 1. Background model
"""

# ╔═╡ 57ee9f4a-fa6a-41e5-a2ce-00107a278ea2
md"""
# The first step in calculating the sensitvity is to have a proper background model, which defines how sensitive our experiment can be!

!!! note "Obtaining the background model"
	To obtain the background model used in this analysis I performed a series of **simulations** of various background **sources** (identified mostly from NEMO-3). 
"""

# ╔═╡ 37966182-3b19-4a91-816f-d576c85dd9ba
md"""
# Simulation setup:
1. Falaise version 5.1.5 
2. Reconstruction pipeline: `MockCalibration -> TKrec -> ChargedParticleTracker -> GammaClusterizer -> SNCuts -> MiModule`
3. Magnetic field **off**
4. 8% FWHM energy resolution
5. Foil geometry: `RealisticFlat`

# Data-cuts (using SNCuts)

- 2 tracks
- 2 foil vertices
- 2 **distinct** associated calo-hits
-  $$E_{sum} \in (300, 3500)$$ keV
- max vertex distance on foil $$r < 50$$ mm
- ToF: $$P_{int} \geq 4%$$ & $$P_{ext} \leq 1%$$

"""

# ╔═╡ 4e6e0813-f5b3-4138-948d-7997bb371cdc
md"""

# Background sources

Simulated 100M events each:
1.  $$2\nu\beta\beta$$  --> $$T^{2\nu}_{1/2} =$$ $(round(SNparams["SeThalf2nu"] / (3600*24*365), sigdigits =2)) yr
2.  $$^{214}$$Bi  --> A = $(round(BkgActivityParams[:Bi214_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  $$^{208}$$Tl  --> A = $(round(BkgActivityParams[:Tl208_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  $$^{40}$$K  --> A = $(round(BkgActivityParams[:K40_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  $$^{234m}$$Pa  --> A = $(round(BkgActivityParams[:Pa234m_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  Radon  --> A = $(round(BkgActivityParams[:Bi214_wire_surface], sigdigits = 3)) $$Bq/m^3$$


"""

# ╔═╡ 83b19179-efaa-4b92-aa71-6cfcea2b26e9
md"""
# Full model
![background model](https://github.com/Shoram444/SNSensitivityEstimate/blob/main/plots/LSM_report/background_model/sumE/background_model_Bfield=Boff_Eres=8perc_mode=sumE_neutron_config=no_neutron_signal=bb0nu_trackAlgo=TKrec.png?raw=true)
"""

# ╔═╡ 49ec8eb7-18ac-4d1e-9b45-401ee2563396
md"""
# Background by source 

| **process**      | **bkg counts** |
|:-----------------|:----------------------|
| 2nubb            | 87663.0±29.0        |
| Bi214            | 8.903±0.022         |
| Radon            | 10.62±0.15          |
| Tl208            | 1.0227±0.0034       |
| K40              | 12155.0±63.0        |
| Pa234m           | 15578.0±39.0        |
| total            | 115417.0±80.0       |


!!! update " "
	- Largest contribution due to $$2\nu\beta\beta$$ - almost 90k events (+10k compared to old tracking)
	  - Good for 2nu physics
      - Bad for 0nu
	- K40 and Pa234m relevant **only** at lower energies 
	- Radon will highly depend on whether or not we have anti-Rn
	- For now, no neutron data

"""

# ╔═╡ b95b52fe-f516-4a41-bdc1-d667ebfbf7bf
md"""
# 2. Methodology
"""

# ╔═╡ 7b318089-8d66-46d9-a829-f74f27b9e565
md"""
# 2.1. Feldman Cousins sensitvity

```math 
T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mathcal{S}(\bar{b})} 
```
"""

# ╔═╡ 48652f8f-274d-4899-b598-1c03b5927961
md"""
# The one you know: 1D approach!

I already gave several talks on this: DocDB#5943 and DocDB#5833

!!! note "The general idea"
	1. Simulate relevant signal process
	2. Chose 1 variable which is tracked (i.e.$$E_{sum}$$, $$E_i$$, $$\phi$$..)
	3. For each possible combination of ROI calculate signal $$\varepsilon$$ and $$\mathcal{S}(\bar{b})$$* as signal-to-background ratio $$r(ROI) = \frac{\varepsilon}{\mathcal{S}(\bar{b})}$$
	4. Find which ROI maximizes $$r(ROI)$$ --> from the $$max(r)$$ calculate sensitvity (simulation)
	5. Once you've measured data --> fit the in ROI to get $$\bar{b}$$ --> calculate sensitivity 
	   * use $$\varepsilon$$ from simulation, $$\bar{b}$$ from data

\* $$\mathcal{S}(\bar{b})$$ is a Fledman Cousins limit calculated at $$\alpha$$ CL for total expected number of background counts in ROI $$\bar{b}$$
"""

# ╔═╡ f5479dd2-e121-4c11-8c63-e31aaeecffca
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

- ROI in this case is $$E_{sum}^l, E_{sum}^u$$, where $$l, u$$ are the lower and upper boundaries of ROI
- Then sensitivity is a function of $$E_{sum}$$ ROI:
```math
T^{1/2}(ROI) = T^{1/2}(E_{sum}^l, E_{sum}^u) = const.\frac{\varepsilon(E_{sum}^l, E_{sum}^u)}{\mathcal{S}(\bar{b}(E_{sum}^l, E_{sum}^u))}
```

"""

# ╔═╡ Cell order:
# ╟─5479bdcf-bb42-4910-b061-befdb1b86b37
# ╟─64580852-4944-45cb-9b9f-b087edb1c7ff
# ╟─32c435b2-ee35-4a65-ba59-b23ee5ba7992
# ╟─432c7c60-37ba-4ade-8935-500652ddc0c3
# ╟─0bd50ff5-d407-4ba1-b4ab-f5fc52ee3c17
# ╟─137994dd-199b-4531-98f0-ba29a44bc4db
# ╟─cea77a5a-92f3-412d-8e50-69f4252eb52a
# ╟─57ee9f4a-fa6a-41e5-a2ce-00107a278ea2
# ╟─37966182-3b19-4a91-816f-d576c85dd9ba
# ╟─4e6e0813-f5b3-4138-948d-7997bb371cdc
# ╟─83b19179-efaa-4b92-aa71-6cfcea2b26e9
# ╟─49ec8eb7-18ac-4d1e-9b45-401ee2563396
# ╟─b95b52fe-f516-4a41-bdc1-d667ebfbf7bf
# ╟─7b318089-8d66-46d9-a829-f74f27b9e565
# ╟─48652f8f-274d-4899-b598-1c03b5927961
# ╠═f5479dd2-e121-4c11-8c63-e31aaeecffca
