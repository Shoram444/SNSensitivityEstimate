### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 96f0e67c-16af-11f0-08a5-ffefdd2c8429
using DrWatson

# ╔═╡ 39a6343a-b95c-4ea4-98e2-83eab6d6c602
@quickactivate "SNSensitivityEstimate"

# ╔═╡ 7d6e5825-c187-4d55-b877-1e1193190ff2
begin
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
end

# ╔═╡ c38eeb3e-11c2-40d0-98cb-89711b799441
using CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, FHist, PrettyTables, StatsBase, ColorSchemes, PlutoUI, CSV

# ╔═╡ 5c201d62-6250-4289-ad8a-bcc519ef570a


# ╔═╡ 33db2852-a205-4a2c-8744-037f5b7a6a80
html"<button onclick='present()'>present</button>"


# ╔═╡ 107241e3-a8b7-43af-b26f-b2b205f7e294
md"""
# Updates on the sensitivity study to BSM physics.
### New and improved multi-dimensional frequentist limits & improved Bayesian models

#### Maros Petro
"""

# ╔═╡ 838099af-cd31-480e-8952-9eda9df9e10a
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

# ╔═╡ 65b97922-2625-43a4-bf9f-622cce900bcf
md"""
# 1. Background model
"""

# ╔═╡ 4011c3d8-d298-40f1-9986-77e991a25152
md"""
# The first step in calculating the sensitvity is to have a proper background model, which defines how sensitive our experiment can be!

!!! note "Obtaining the background model"
	To obtain the background model used in this analysis I performed a series of **simulations** of various background **sources** (identified mostly from NEMO-3). 
"""

# ╔═╡ 781740e1-5f0c-4b99-834c-1ac6a75118b7
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

# ╔═╡ 4db07780-248f-4658-82a5-507ea65edb0b
md"""

# Background sources

Simulated 100M events each:
1.  $$2\nu\beta\beta$$  --> $$T^{2\nu}_{1/2} =$$ $(round(SNparams["SeThalf2nu"] / (3600*24*365), sigdigits =2)) yr
2.  $$^{214}$$Bi  --> A = $(round(BkgActivityParams[:Bi214_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  $$^{208}$$Tl  --> A = $(round(BkgActivityParams[:Tl208_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  $$^{40}$$K  --> A = $(round(BkgActivityParams[:K40_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  $$^{234m}$$Pa  --> A = $(round(BkgActivityParams[:Pa234m_foil_bulk], sigdigits = 3)) $$Bq/kg$$
2.  Radon  --> A = $(round(BkgActivityParams[:Bi214_wire_surface], sigdigits = 3)) $$Bq/m^3$$

🤔 Neutron data will be added later. I received the files from Sam but haven't yet been able to add them to my framework. 🤭

"""

# ╔═╡ 16fbb9bb-d86b-433b-91fa-a12f8008eaa9
md"""
# Background Model
"""

# ╔═╡ 5099a914-3f46-4e86-9570-bc13f5d3ca80
md"""
Emin = $(@bind Emin PlutoUI.Slider(300:100:3500, default=0; show_value=true)) keV

Emax = $(@bind Emax PlutoUI.Slider(0:100:3500, default=3500; show_value=:true)) keV

Thalf = $(@bind Thalf PlutoUI.Slider((1e3:-200:1)*1e24, default=1e27; show_value=:true)) yr
"""

# ╔═╡ 1f4ae42b-61be-44ed-bb82-60e7cfb582ba
begin
	files_directory = "fal5_8perc_Boff_TKrec_evis_bcu_J38"

	data_processes = load_data_processes(
	    files_directory, 
	    "sumE",
	    fwhm = 0.0
	)
	signal = get_process("bb0nu_foil_bulk", data_processes)

    background = [
        get_process("bb_foil_bulk", data_processes),
        get_process("Bi214_foil_bulk", data_processes),
        get_process("Bi214_wire_surface", data_processes),
        get_process("Tl208_foil_bulk", data_processes),
        get_process("K40_foil_bulk", data_processes),
        get_process("Pa234m_foil_bulk", data_processes),
    ]
	set_nTotalSim!( signal ,1e8)
	set_nTotalSim!( background[1], 0.99e8 )
    set_nTotalSim!( background[2], 1e8 )
    set_nTotalSim!( background[3], 1e8 )
    set_nTotalSim!( background[4], 1e8 )
    set_nTotalSim!( background[5], 1e8 )
    set_nTotalSim!( background[6], 1e8 )
	set_signal!(background[1], false)
	nothing
end

# ╔═╡ de64766d-fda7-4a0e-bc7c-80bb90554bac
begin
	# Total background model:
    bkg_hists = [restrict(get_bkg_counts_1D(b), Emin, Emax) for b in background]
    sig_hist = get_bkg_counts_1D(signal)


    with_theme(theme_latexfonts()) do
        f = Figure()
        ax = Axis(
            f[1,1], 
            # xlabel = analysisDict[:mode], 
            xlabel = L"$E_1 + E_2$ (keV)", 
            ylabel = L"counts / $17.5$kg.yr exposure / $100$ keV", 
            yscale = log10, 
            limits = (Emin, Emax, 1e-4, 1e6),
            # limits = (2500, 3500, 0, 3),
            title = "Total background model\nsummed 2-electron energy"
        )
        
        colors = ColorSchemes.tol_vibrant
        labels = [b.isotopeName for b in background]
        st = hist!(ax, sum(bkg_hists), label =labels[1],color=colors[1], strokewidth = 1, strokecolor = :black)
        errorbars!(ax, sum(bkg_hists), color = :black, whiskerwidth = 7)
        
        for i=2:length(bkg_hists)
            hist!(ax, sum(bkg_hists[i:end]), label=labels[i], color=colors[i], strokewidth = 1, strokecolor = :black)
            
        end
        lines!(ax, midpoints(binedges(sig_hist)), bincounts(sig_hist) .* (Thalf/1e26), label = signal.isotopeName, color = :black, linestyle = :dash, linewidth = 2.5)

        ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
        ax.xticks = 0:500:3500
        Legend(f[2,1], ax, orientation=:horizontal, fontsize=8, nbanks = 3)
        f
    end
end

# ╔═╡ 28a3b440-03b7-46fc-a08c-46084b41ae17
bkg_hists[1] |> integral

# ╔═╡ bad82399-f4f0-401c-acc7-b7f67155f27f
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

# ╔═╡ 4a2ca349-25ff-4a16-a4cf-075767700240
md"""
# 2. Methodology
"""

# ╔═╡ eb22cfcd-a08e-4e57-9f02-113612804fc1
md"""
# 2.1. Feldman Cousins sensitvity

```math 
T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mathcal{S}(\bar{b})} 
```
"""

# ╔═╡ aae87c04-2352-4db5-9018-f8d0fa75e4b9
md"""
# The one you know: 1D approach!

I already gave several talks on this: DocDB#5943 and DocDB#5833

!!! note "The general idea"
	1. Simulate relevant signal process
	2. Chose 1 variable which is tracked (i.e.$$E_{sum}$$, $$E_i$$, $$\phi$$..)
	3. For each possible combination of ROI calculate signal $$\varepsilon$$ and $$\mathcal{S}(\bar{b})$$* as signal-to-background ratio $$r(ROI) = \frac{\varepsilon}{\mathcal{S}(\bar{b})}$$
	4. Find which ROI maximizes $$r(ROI)$$ --> from the $$max(r)$$ calculate sensitvity (simulation)
	5. Once you've measured data --> fit in ROI to get $$\bar{b}$$ --> calculate sensitivity 
	   * use $$\varepsilon$$ from simulation, $$\bar{b}$$ from data

\* $$\mathcal{S}(\bar{b})$$ is a Fledman Cousins limit calculated at $$\alpha$$ CL for total expected number of background counts in ROI $$\bar{b}$$
"""

# ╔═╡ c8539419-4919-48a9-ad41-0dcba22318e8
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

- ROI in this case is $$E_{sum}^l, E_{sum}^u$$, where $$l, u$$ are the lower and upper boundaries of ROI
- Then sensitivity is a function of $$E_{sum}$$ ROI:
```math
T^{1/2}(ROI) = T^{1/2}(E_{sum}^l, E_{sum}^u) = const.\frac{\varepsilon(E_{sum}^l, E_{sum}^u)}{\mathcal{S}(\bar{b}(E_{sum}^l, E_{sum}^u))}
```
- We calculate $$T^{1/2}$$ for each combination and create a 2D map:



"""

# ╔═╡ c519cbff-afb1-4ac7-9bfe-1e07cf78979b
begin
	global α = 1.64485362695147
	t12MapESum = get_tHalf_map(SNparams, α, signal, background...;)
    best_t12ESum = get_max_bin(t12MapESum)

    # If you want additional info, like background counts in the ROI, use: `get_bkg_counts_ROI`.
    expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)

    # To get the signal efficiency at the ROI, use: `lookup(signal, best_t12ESum)`
    effbb = lookup(signal, best_t12ESum)
	ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)
	lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
          b  = $(round(expBkgESum, digits = 2)) 
          T12 ≥  $(ThalfbbESum) yr 
          ε = $(round(effbb, digits = 2)*100)%"

    let 
        f = Figure(size=(600, 400))
        a = Axis(f[1,1], xlabel = "min_ROI [keV]", ylabel = "max_ROI [keV]")
        p = plot!(a, t12MapESum)
        text!(a, 2000, 500, text=lbl)
        Colorbar(f[1,2], p, label="sensitivity [yr]", scale=log10)
        f
    end
end

# ╔═╡ d7354883-f530-4bfe-955d-ca37263d0c7b
md"""
- best ROI is $$2700 - 3100$$ keV
- with $$\bar{b}=1.04$$ and $$\varepsilon = 15%$$ we get sensitivity of $$T^{1/2}\geq 4.14\cdot 10^{24}$$ yr

!!! danger "What if we look at another variable? Another signal process?"
	- For different signal processes (signal shapes!) we can be more sensitive in different channels: i.e. $$\phi$$ **should** be a better channel for RH spectra (will get to this later)
	- Do we really need to chose only one? What if we can maximize **all of them at once!**

"""

# ╔═╡ 87c50992-3347-47e8-8cb4-3f0dcee41ebf
md"""
# N-Dim ROI Search: the next step!

!!! note "The general idea"
	1. Simulate relevant signal process (again...)
	2. Chose n variables, i.e. $$E_{sum}$$, $$\phi$$, $$r$$
	   - That means we have 6 parameters (2 for each ROI bound)
	   - n can be essentially any number, but computational requirements grow fast (Can be used to study optimal vertex-distance, $$E_{max}$$, $$E_{min}$$, 2D angle etc.)
	3. For each possible combination of ROI calculate signal $$\varepsilon$$ and $$\mathcal{S}(\bar{b})$$* as signal-to-background ratio $$r(ROI) = \frac{\varepsilon}{\mathcal{S}(\bar{b})}$$ 
	   - 💀**This takes way too long and grows way too fast, when we add new dimensions** 💀
	   - 🍰**We instead more advanced optimization algorithms that cut down the computing time!** 🍰
	4. Find which ROI maximizes $$r(ROI)$$ --> from the $$max(r)$$ calculate sensitvity (simulation)
	5. Once you've measured data --> fit in ROI to get $$\bar{b}$$ --> calculate sensitivity 
	   * use $$\varepsilon$$ from simulation, $$\bar{b}$$ from data

"""

# ╔═╡ 6fa76137-73fc-47ea-927b-0ffed51149f7
md"""
# Example for $$0\nu\beta\beta$$ N-Dim:

- Choose n = 3, with variables $$\varphi, E_{sum}, r$$
  - ROI consists of 6 parameters: ($$E_{sum}^l, E_{sum}^u$$), ($$\varphi^l, \varphi^u$$) and ($$r^l, r^u$$), where $$l, u$$ are the lower and upper boundaries of ROI
   - let's collectively call them $$\Theta$$
- Choose optimization algorithm 
  - there's many libraries written for this purpose, the goal is to **choose one that works well for the particular problem**
  - I tried multiple, the best results seem to be found with **Evolutionary Centers Algorithm** [ECA](https://link.springer.com/chapter/10.1007/978-981-13-0860-4_6)¹
- Then sensitivity is a function of $$\Theta$$:
```math
T^{1/2}(\Theta) = const.\frac{\varepsilon(\Theta)}{\mathcal{S}(\bar{b}(\Theta))}
```


¹https://doi.org/10.1007/978-981-13-0860-4_6
"""

# ╔═╡ cdcc36a8-888b-42e4-8721-18e183a3ae5e
md"""
# Example for $$0\nu\beta\beta$$ N-Dim:

!!! update "Basics of ECA approach"
	1. The algorithm starts with N randomly placed initial guesses $$\Theta_i$$
	2. Calculates the values of $$T^{1/2}(\Theta_i)$$ for each
	3. Identifies "best performing" $$\Theta_i$$
	4. Creates new set of N around the best performers
	5. Repeat until convergence (or failure...)
"""

# ╔═╡ 99c39ae6-9196-4b28-be41-fc885ccf7265
md"""
# Example for $$0\nu\beta\beta$$ N-Dim:

Plots of E, phi, r scatters
"""

# ╔═╡ 724e4ff0-6dea-4a73-8708-0c886c98a8db
md"""
# Example for $$0\nu\beta\beta$$ N-Dim:

Results
$$\Theta =$$

|variable | ROI |
|:-------:|:----:|
|$$E_{sum}$$ |(2700, 3200) keV|
|$$E_{max}$$| (1000, 3000) keV|
|$$\varphi$$| (0, 180)$$\degree$$|
|$$r$$ | (0, 60) mm |

And 

|variable | value |
|:-------:|:----:|
|$$\varepsilon$$ |0.15|
|$$\bar{b}$$| 1.04 ± 0.02|
|$$\mathcal{S}$$| 3.3|
|$$T^{1/2}$$ | $$\geq 4.04 yr$$|

"""

# ╔═╡ 8680f3a5-6ed9-40de-b8a3-c2e34b4507e0
md"""
# 2.2 Bayesian sensitivity

```math 
T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_S} 
```
where $$\mu_S$$ is the posterior 90% CI on the signal!
"""

# ╔═╡ e6bbd64e-a2eb-462f-837a-079ef07f3c43
md"""
# 2.2 Bayesian sensitivity

!!! note "The general idea"
	1. Simulate relevant signal and background processes
       - form **truth** pdfs on signal and background processes
	2. Chose 1 variable which is tracked (i.e.$$E_{sum}$$, $$E_i$$, $$\phi$$..)
	3. Create appropriate model which describes the data
	   - One of the model's parameters is $n_S$ -> the expected number of signal events
	4. Generate a pseudo-data from the simulated pdfs 
	5. Extract **pdf** of the signal counts **given data**: $$pdf(n_S|data)$$ 
	6. Extract the 90% CI from $$pdf(n_S|data)$$ --> Calculate $$T^{1/2}$$
	7. Repeat steps 4.-6. N times to mitigate statistical fluctuations and obtain **median sensitivity**
	8. Profit 🍾🍾🍾
	
More detailed description in: DocDB#5943 
"""

# ╔═╡ 085a4b0a-1fdc-44d8-a8ce-7158ac54256f
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:
!!! danger "What was...bad!"
	In DocDB#5943 I showed **general model** 
	**Exponential background + gaussian signal** 
	```math
		\mathcal{L}(data|n_S, n_B, \lambda) = \frac{1}{n_S + n_B}\prod_i^{N^{obs}}\left( n_B\lambda e^{-\lambda E_i} + n_S\frac{1}{\sigma\sqrt{2\pi}}e^{-(\frac{E_i-\mu}{\sigma})^2} \right)
	```
	Where $$\mathcal{L}(data|n_S, n_B, \lambda)$$ is the likelihhod of data **given** parameters, $$n_S, n_B, \lambda$$ are expected number of signal, background counts and decay constant. $$N_{obs}$$ are the observed events.

	This description works for most experiments **that are not SuperNEMO!!**

	- Our ROI is very wide:
	   - cannot be described well with and exponential, or flat (as is the case for CUPID and Gerda)
	- Our Signal is not gaussian:
	   - due to energy losses in gas, the signal is more Landau shaped!

"""

# ╔═╡ 187dfcf3-48d8-4ee9-8d33-8b3bac063b37
md"""
!!! note "What is...better! 🥁"
	**Detailed** likelihood model with each process having its own pdf!
	```math
			\mathcal{L}(data|\vec{\Theta}) = \frac{1}{\sum\Theta_i}\prod_i^{N^{obs}}\left( \Theta_1 pdf(sig, E_i) + \sum_{j=2}^{n_{pars}}\Theta_j pdf(bkg_j, E_i) \right)
	```
	Where $$pdf(sig, E_i)$$ and $$pdf(bkg_j, E_i)$$ are the signal and background normalized spectra evaluated at $$E_i$$. 

	- This way we have better description of the shapes of the spectra (if we consider simulation to be correct...)
"""

# ╔═╡ 35e1a3ac-3f0d-482e-9e19-fdc00a0a9004


# ╔═╡ 87d20c02-cd0a-43a1-a437-016086bf76fd


# ╔═╡ bf7e38c5-d254-4392-b526-081784269fd8


# ╔═╡ e632053f-fbc5-41ab-a9d5-c3a826919402


# ╔═╡ Cell order:
# ╠═5c201d62-6250-4289-ad8a-bcc519ef570a
# ╠═33db2852-a205-4a2c-8744-037f5b7a6a80
# ╠═96f0e67c-16af-11f0-08a5-ffefdd2c8429
# ╠═39a6343a-b95c-4ea4-98e2-83eab6d6c602
# ╠═c38eeb3e-11c2-40d0-98cb-89711b799441
# ╠═7d6e5825-c187-4d55-b877-1e1193190ff2
# ╠═107241e3-a8b7-43af-b26f-b2b205f7e294
# ╠═838099af-cd31-480e-8952-9eda9df9e10a
# ╠═65b97922-2625-43a4-bf9f-622cce900bcf
# ╠═4011c3d8-d298-40f1-9986-77e991a25152
# ╠═781740e1-5f0c-4b99-834c-1ac6a75118b7
# ╠═4db07780-248f-4658-82a5-507ea65edb0b
# ╠═16fbb9bb-d86b-433b-91fa-a12f8008eaa9
# ╠═5099a914-3f46-4e86-9570-bc13f5d3ca80
# ╠═de64766d-fda7-4a0e-bc7c-80bb90554bac
# ╠═1f4ae42b-61be-44ed-bb82-60e7cfb582ba
# ╠═28a3b440-03b7-46fc-a08c-46084b41ae17
# ╟─bad82399-f4f0-401c-acc7-b7f67155f27f
# ╟─4a2ca349-25ff-4a16-a4cf-075767700240
# ╟─eb22cfcd-a08e-4e57-9f02-113612804fc1
# ╟─aae87c04-2352-4db5-9018-f8d0fa75e4b9
# ╟─c8539419-4919-48a9-ad41-0dcba22318e8
# ╟─c519cbff-afb1-4ac7-9bfe-1e07cf78979b
# ╟─d7354883-f530-4bfe-955d-ca37263d0c7b
# ╟─87c50992-3347-47e8-8cb4-3f0dcee41ebf
# ╟─6fa76137-73fc-47ea-927b-0ffed51149f7
# ╠═cdcc36a8-888b-42e4-8721-18e183a3ae5e
# ╠═99c39ae6-9196-4b28-be41-fc885ccf7265
# ╠═724e4ff0-6dea-4a73-8708-0c886c98a8db
# ╠═8680f3a5-6ed9-40de-b8a3-c2e34b4507e0
# ╠═e6bbd64e-a2eb-462f-837a-079ef07f3c43
# ╠═085a4b0a-1fdc-44d8-a8ce-7158ac54256f
# ╠═187dfcf3-48d8-4ee9-8d33-8b3bac063b37
# ╠═35e1a3ac-3f0d-482e-9e19-fdc00a0a9004
# ╠═87d20c02-cd0a-43a1-a437-016086bf76fd
# ╠═bf7e38c5-d254-4392-b526-081784269fd8
# ╠═e632053f-fbc5-41ab-a9d5-c3a826919402
