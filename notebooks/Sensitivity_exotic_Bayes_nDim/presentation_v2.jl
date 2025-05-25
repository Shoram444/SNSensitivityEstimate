### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
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

# ╔═╡ 992a6366-86e0-4d00-ab9c-ac29190ff693
begin
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
end

# ╔═╡ 6d3bd21c-a856-407e-b3e7-d98d8cad0c9b
using Distributions

# ╔═╡ c38eeb3e-11c2-40d0-98cb-89711b799441
using CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, FHist, PrettyTables, StatsBase, ColorSchemes, PlutoUI, CSV

# ╔═╡ 7d6e5825-c187-4d55-b877-1e1193190ff2
include(srcdir("params/Params.jl"))

# ╔═╡ 1a041346-47c7-41a4-a30c-5710a5dd8786
html"""
<style>
	pluto-output { font-size : 24px;}

main {
	margin: 0 auto;
	max-width: 2000px;
	padding-left: max(160px, 5%);
	padding-right: max(160px, 5%);
}

</style>
"""

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
### 2.  Methodology
   1. #### Feldman Cousins sensitivity
   2. #### Bayesian sensitivity
### 3. Results
   1. #### $$0\nu\beta\beta$$
   2. #### $$0\nu\beta\beta\chi^0$$ 
   3. #### $$0\nu\beta\beta\chi^0\chi^0$$
   4. #### $$\nu_R\nu_L\beta\beta$$
   
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
2. Reconstruction pipeline: `MockCalibration -> ` **Cimrman** `-> ChargedParticleTracker -> GammaClusterizer -> SNCuts -> MiModule`
3. Magnetic field **off**
4. 8% FWHM energy resolution
5. Foil geometry: `RealisticFlat`

# Data-cuts (using SNCuts)

- 2 tracks
- 2 foil vertices
- 2 **distinct** associated calo-hits
-  $$E_{sum} \in (300, 3500)$$ keV
- max vertex distance on foil $$r < 50$$ mm
- ToF: $$P_{int} \geq 4\%$$ & $$P_{ext} \leq 1\%$$

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

Thalf = $(@bind Thalf PlutoUI.Slider((1:10:1000)*1e22, default=1e24; show_value=:true)) yr
"""

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

# ╔═╡ c12114cf-0b28-4579-8bbb-147e8acc6f33
html"""
<style>
pluto-output table > tbody td {
    font-size: 20px;
}

</style>
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

# ╔═╡ 6edb8cab-331c-403c-8c93-b8af84f38989
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

Boundaries of ROI affect signal + background differently!

That we are interested in is:
- signal efficiency $$\varepsilon(ROI)$$
- expected background counts $$\bar{b}(ROI)$$
"""

# ╔═╡ 1e127dad-951d-4873-b711-fb8a23f3d712
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

For each combination of ROI ($$E_{sum}^l, E_{sum}^u$$) we get a different value of $$T^{1/2}$$ --> pick maximum
"""

# ╔═╡ 755fcdeb-d7a8-4d41-bb7d-2ec9b5e3065b
md"""
- best ROI is $$2700 - 3100$$ keV
- with $$\bar{b}=1.04$$ and $$\varepsilon = 15%$$ we get sensitivity of $$T^{1/2}\geq 4.14\cdot 10^{24}$$ yr

"""

# ╔═╡ d7354883-f530-4bfe-955d-ca37263d0c7b
md"""
#  
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
	   - 🍰**We instead use more advanced optimization algorithms that cut down the computing time!** 🍰
	4. Find which ROI maximizes $$r(ROI)$$ --> from the $$max(r)$$ calculate sensitvity (simulation)
	5. Once you've measured data --> fit in ROI to get $$\bar{b}$$ --> calculate sensitivity 
	   * use $$\varepsilon$$ from simulation, $$\bar{b}$$ from data

"""

# ╔═╡ 9ae6e7cf-fa54-4d17-834f-a11f7201bdff
md"""
# Extending to multiple dimensions:
Let's take the example methodology above, but use 2 variables for best ROI: $$E_{sum}$$ and $$\varphi$$

**Again, we look at the spectrum for signal and background (now a 2D histogram):**
"""

# ╔═╡ 112a6be9-e6ac-4ee2-b7b4-82fc7d5d7b88
md"""
 $$E^l_{sum}$$ = $(@bind E_l PlutoUI.Slider(0:100:3500, default =2500, show_value = true)),   
 $$E^u_{sum}$$ = $(@bind E_u PlutoUI.Slider(0:100:3500, default =3300, show_value = true)),
 $$\varphi^l$$ = $(@bind p_l PlutoUI.Slider(0:5:180, default =10, show_value = true)),
 $$\varphi^u$$ = $(@bind p_u PlutoUI.Slider(0:5:180, default =160, show_value = true))
"""

# ╔═╡ 8bc26abf-c593-44d8-b1d1-5f2db8d07705
md"""
We can see that taking a 2D ROI (basically a data-cut), we can keep a lot of signal while reducing the amount of background.
In this example however, it is not exactly clear which is the **best** choice!


!!! note ""
	This approach can be extended to many dimensions!
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

# ╔═╡ 724e4ff0-6dea-4a73-8708-0c886c98a8db
md"""
# Example for $$0\nu\beta\beta$$ N-Dim:

```math
\textrm{Best} ~~\Theta = 
```

|variable | ROI |
|:-------:|:----:|
|$$E_{sum}$$ |(2710, 3350) keV|
|$$\varphi$$| (10, 180)$$\degree$$|
|$$r$$ | (0, 50) mm |

```math
\textrm{Resulting sensitivity}
```
|variable | value |
|:-------:|:----:|
|$$\varepsilon$$ |0.145|
|$$\bar{b}$$| 0.73|
|$$T^{1/2}$$ | $$\geq 4.22 yr$$|

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

# ╔═╡ e768fb37-96d1-4181-9e36-1aa5af23db1f
md"""
# Example, generic Bayesian inference

Let's say we have some Gaussian signal and Exponential background and have measured some data. 

!!! note "We have:"
	*  $$\mu_{sig}$$ mean of signal Gaussian
	*  $$\Theta_{sig}$$ proportion of signal in data
	*  $$\lambda_{bkg}$$ decay rate of background
	*  $$\Theta_{bkg}$$ proportion of background in data


!!! update ""
	* For simplicity, let's say we know $$\mu_{sig}$$ and $$\lambda_{bkg}$$. 
	* We just want to find: $$\Theta_{sig}$$, $$\Theta_{bkg}$$.

\* let's ignore $$\sigma$$ for now..
"""

# ╔═╡ 8c041449-ec1c-4bdf-ac7b-0ef59ba777a4
md"""
# Example, generic Bayesian inference

Our likelihood:
```math
		\mathcal{L}(data|\Theta_{sig}, \Theta_{bkg}) = \prod_i^{N^{obs}}\left( \Theta_{bkg}\lambda_{bkg} e^{-\lambda_{bkg} E_i} + \Theta_{sig}\frac{1}{\sigma\sqrt{2\pi}}e^{-(\frac{E_i-\mu_{sig}}{\sigma})^2} \right)
```

Example experiment with n data-points:
"""

# ╔═╡ 1910efa8-1530-4139-8287-d9d05c3e69ae
md"""
 n = $(@bind n PlutoUI.Slider(5:100:3000, default = 100, show_value = true))
"""


# ╔═╡ e83967c4-4bf9-4715-8aca-d022209c77d5
md"""
Parameters of the experiment:
"""

# ╔═╡ 518693c6-694b-45eb-a08b-c3a449a0d0bf
md"""
 $$\mu_{sig}$$ = $(@bind μ_sig PlutoUI.Slider(1.0:1.0:10.0, default = 5.0, show_value = true)) 
"""

# ╔═╡ ea51b339-1d49-4f04-a8fe-744f9e855026
md"""
 $$\Theta_{sig}$$ = $(@bind n_sig PlutoUI.Slider(0.001:0.1:1.0, default = 0.1, show_value = true))
"""

# ╔═╡ ed6b7268-412d-45a9-895c-7248b0686efe
md"""
 $$\lambda_{bkg}$$ = $(@bind lambda_bkg PlutoUI.Slider(0.5:0.1:3.0, default = 0.8, show_value = true))
"""

# ╔═╡ c57ec8c2-68d1-47b9-b7e7-0a415ce85031
md"""
 $$\Theta_{bkg}$$ = 1 - $$\Theta_{sig}$$ = $(n_bkg = 1.0-round(n_sig, sigdigits= 1))
"""

# ╔═╡ 5b254b0e-6b4b-46be-bed9-8df0efa3de44
begin
	
	signal_data = rand(Normal(μ_sig, 1.0), Int(round(n*n_sig)))
	background_data = rand(Exponential(lambda_bkg), Int(round(n*n_bkg)))

	h1_sig = Hist1D(signal_data; binedges = 0:0.2:10.0)
	h1_bkg = Hist1D(background_data; binedges = 0:0.2:10.0)

	true_ps = n*n_sig # true proportion of signal
	true_pb = n*n_bkg
	signal_dist = Normal(μ_sig, 1.0)
	background_dist = Exponential(lambda_bkg)
	
	data = vcat(signal_data, background_data)

	f_ex = Figure()
	ax_ex = Axis(f_ex[1,1], xlabel = "arbitrary energy", ylabel = "counts")
	xs = binedges(h1_sig) |> midpoints |> collect
	p = lines!(ax_ex, xs, x-> pdf(Normal(μ_sig, 1.0), x) .* n .* 0.2 , label = "signal", color = :blue, linewidth = 4)
	lines!(ax_ex,  xs, x-> pdf(Exponential(lambda_bkg), x) .* n .* 0.2, label = "background", color= :red, linewidth = 4)
	stairs!(ax_ex, h1_sig + h1_bkg, label = "data", linewidth = 4)
	axislegend(ax_ex)
	f_ex
end

# ╔═╡ 9043a542-5102-488e-a3ab-0f377648145a
md"""
For a Bayesian inference we must chose the so-called **prior** (describes our knowledge of the paramaters)

```math
	pdf(\Theta) \sim \textrm{some distribution}
```

What we really want as the result is the so-called **posterior**: 
**The distribution of the parameter given data!**

Since our parameter of interest is $$\Theta_{sig}$$:
```math
	pdf(\Theta_{sig}|data) \propto \mathcal{L}(data|\Theta_{sig}, \Theta_{bkg}) \cdot pdf(\Theta_{sig})
```

"""

# ╔═╡ e5664e79-6d50-43dc-963a-d2207c48b9a8
md"""
p1 = $(@bind prior_par1 PlutoUI.Slider(1:1:10.0, default = 1, show_value = true))
"""

# ╔═╡ cc0d49b8-3d6b-4a2c-8300-0aea834df028
md"""
p2 = $(@bind prior_par2 PlutoUI.Slider(1:1:5.0, default = 1, show_value = true))
"""

# ╔═╡ 25f1adc4-a6f4-4a10-9fe1-b0d1817f0ece
begin
	# --- Grid for p_s ---
	ps_grid = range(0.0, 1.0; length=500)
	log_posterior = zeros(length(ps_grid))
	log_likelihood_vector = zeros(length(ps_grid))

	prior_signal = Beta(prior_par1,prior_par2)
	# --- Evaluate likelihood over grid ---
	for (i, ps) in enumerate(ps_grid)
	    pb = 1 - ps
	    log_likelihood = sum(log.(ps * pdf.(signal_dist, data) .+ pb * pdf.(background_dist, data)))
		log_prior = logpdf(prior_signal, ps)
		log_likelihood_vector[i] = log_likelihood
	    log_posterior[i] = log_likelihood + log_prior # Flat prior assumed
	end

	
	
	# --- Normalize to get posterior ---
	posterior = exp.(log_posterior .- maximum(log_posterior))  # for numerical stability
	posterior ./= sum(posterior)

	fit_f(x, p) = p[1]*pdf(Normal(μ_sig, 1.0), x) + p[2]*pdf(Exponential(lambda_bkg), x)

	p_best = [argmax(posterior)/length(posterior), 1-argmax(posterior)/length(posterior)]
	
	let 
	f = Figure(size = (1000, 700))
	ax = Axis(f[1,1:2], xlabel = "arbitrary energy", ylabel = "counts")
	xs = binedges(h1_sig) |> midpoints |> collect
	stairs!(ax, h1_sig + h1_bkg, label = "data", linewidth = 4)
	p = lines!(ax, xs, x-> pdf(Normal(μ_sig, 1.0), x) .* n .* 0.2, label = "signal", color = :blue, linewidth = 4)
	lines!(ax,  xs, x-> pdf(Exponential(lambda_bkg), x) .* n .* 0.2, label = "background", color= :red, linewidth = 4)
	lines!(ax,  xs, x-> fit_f(x, p_best) .* n .* 0.2, label = "fit", color= :purple, linewidth = 4)
	axislegend(ax)


	a =Axis(f[3,1:2], xlabel = L"\Theta", ylabel = L"pdf($\Theta_i$| data)", title= "posterior distribution")
	p = lines!(a, ps_grid, posterior, label = L"\Theta_{sig}")
	lines!(a, 1.0 .- ps_grid, posterior, label = L"\Theta_{bkg}")
	axislegend(a)

	a_prior = Axis(f[2,2], title = "prior signal")
	lines!(a_prior, ps_grid, x-> pdf(prior_signal, x))
	
	a_ll = Axis(f[2,1], title = "log likelihood", limits = (nothing, nothing, maximum(log_likelihood_vector)*1.5,maximum(log_likelihood_vector)*0.9))

	lines!(a_ll, ps_grid, log_likelihood_vector, label = "log_likelihood")
	scatter!(a_ll, argmax(log_likelihood_vector)/length(ps_grid), maximum(log_likelihood_vector), label = "mle value = $(argmax(log_likelihood_vector)/length(ps_grid))", color = :black)
	axislegend(a_ll)
		
	f
	end
end

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
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:
!!! note "What is...better! 🥁"
	**Detailed** likelihood model with each process having its own pdf!
	```math
			\mathcal{L}(data|\vec{\Theta}) = \frac{1}{\sum\Theta_i}\prod_i^{N^{obs}}\left( \Theta_1 pdf(sig, E_i) + \sum_{j=2}^{n_{pars}}\Theta_j pdf(bkg_j, E_i) \right)
	```
	Where $$pdf(sig, E_i)$$ and $$pdf(bkg_j, E_i)$$ are the signal and background normalized spectra evaluated at $$E_i$$. 

	- This way we have better description of the shapes of the spectra (if we consider simulation to be correct...)

!!! update "Priors and posetrior"
	**Priors**
	  - the parameters are set up as a "proportion of the total spectrum", so each process contributes some percentage of the total spectrum: $$\Sigma_i \Theta_i = 1$$
	  - signal $$p(\Theta_1 | data) \sim Uniform(0, 10^{-4})$$ 
	    -  $$10^{-4}$$ represents roughly $$T^{1/2}\geq 10^{24} yr$$ for $$17.5 kg.y$$
	  - backgrounds $$\Theta_i \sim Uniform(0,1)$$

!!! danger ""
	The only condition is that $$\sum \Theta_i = 1.0$$. 


"""

# ╔═╡ 8a0763e0-bf30-49df-930c-d28ded622ded
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

!!! danger "Posterior"
	The posterior distributions represent the **possible distributions of the parameter** given data. That means, that for each parameter (proportion of the spectrum), we get a pdf. 
	For example, if $$0\nu\beta\beta$$ represented say 30% of our measured data, we'd get a pdf with peak around 30%. If there's no signal present, the peak would be expected at 0. The shape of the pdf is how we calculate confidence interval.
"""

# ╔═╡ 559c3769-c956-47e7-b5ea-07b866333104
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

Results:
For a single Bayesian inference we get posterior distribution of $$\Theta_{sig}$$:


"""

# ╔═╡ 4e8c3fa4-b9b4-4656-a61e-2756ce10f3f1
html"""
<img src="https://github.com/Shoram444/SNSensitivityEstimate/blob/main/notebooks/Sensitivity_exotic_Bayes_nDim/samples_As.png?raw=true"
	width="800"
>
"""

# ╔═╡ bad14d75-c550-4221-884c-d57c787f3039
md"""

This means, the most probable value for $$\Theta_{sig}$$ is very close to zero, with 90% CI around $$3\times 10^{-5}$$, which corresponds to 3.5 events. This is the 90% CI on the number of signal events that should be plugged into the sensitivity equation.
"""

# ╔═╡ b2a28a93-2524-43e3-a9fa-6ef30ac213af
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$:

!!! danger "Danger!"
	The Bayesian inference method shown can fluctuate a lot! Both due to the statistical nature of the Bayesian MCMC sampling and due to how the data fluctuates (this is true for Frequentist as well). To limit the effects, we create **n pseudo experiments** to gather a large enough statistics for drawing conclusions. 
"""

# ╔═╡ 8e5fa765-c8a9-4419-9e4f-c58a213be563
md"""
# Example for $$0\nu\beta\beta$$ and $$E_{sum}$$

"""

# ╔═╡ 917eca81-62a3-49e5-a9d1-47d2b6f9c129
html"""
<img src="https://github.com/Shoram444/SNSensitivityEstimate/blob/main/scripts/0nu/Bayes_hist_models/result_bb0nu_foil_bulk_sumE.png?raw=true"
width = 700>
"""

# ╔═╡ a1d02f00-354a-4ebb-b3f0-388c2a69c210
md"""

I ran this algorithm in 20 parallel jobs on CC-LYON cluster for 2 days each.
The following steps were within each job:
1. For each iteration **generate** pseudo-data based on the background model
2. Perform bayesian inference and extract **Sensitivity**
3. Repeat for 2 days
4. Save all sensitivities into a file

The resulting median sensitivity for $$0\nu\beta\beta$$ is:
```math
	T^{0\nu}_{1/2} \geq 4.66 \times 10^24 yr
```
"""

# ╔═╡ 7da00dac-d0a8-4b77-973e-9e46e5fa964d
md"""
# Results
Numbers in parenthasis are with **manually** added "total neutrons in given ROI" (without performing a dedicated ND/Bayes analysis)


|signal|1D|ND|Bayes|
|:----:|:-:|:-:|:-:|
|$$0\nu\beta\beta$$|$$4.14 \times 10^{24} y$$|$$4.22 (3.20) \times 10^{24}y$$|$$4.66 \times 10^{24} y$$|
|$$0\nu\beta\beta\chi^0$$|$$1.45 \times 10^{23} y$$| $$1.48 (1.45) \times 10^{23} y$$ | $$2.38 \times 10^{23} y$$ |
|$$0\nu\beta\beta\chi^0\chi^0$$|$$2.31 \times 10^{22} y$$|$$2.32(2.32) \times 10^{22} y$$  |$$1.43 \times 10^{22} y$$|
|$$\nu_L\nu_R\beta\beta$$| $$1.30 \times 10^{22} y^*$$| $$1.30 \times 10^{22} y$$ | $$1.09 \times 10^{21} y^{**}$$|


 \* for angular distribution

 \*\* for angular distribution, need to investigate why this is so low.  
"""

# ╔═╡ eb626008-caed-4c21-b563-1f0db822355f
md"""
# Results: $$0\nu\beta\beta\chi^0$$
ROI: $$\varphi \in (5, 175)^{\circ}, E_{sum} \in (2500, 3100)keV, r \in (0, 50)mm$$

 $$T^{1/2} \geq 1.48 \times 10^{23}$$ y; $$\varepsilon = 0.017$$; $$\bar{b} = 25.944 ~(+1.2)$$
 
"""

# ╔═╡ e2e4a637-dc3e-4cd5-a156-9cd0ebe56495
md"""
# Results: $$0\nu\beta\beta\chi^0\chi^0$$
ROI: $$\varphi \in (0, 180)^{\circ}, E_{sum} \in (1200, 3100)keV, r \in (0, 50)mm$$

 $$T^{1/2} \geq 2.32 \times 10^{22}$$ y; $$\varepsilon = 0.09$$; $$\bar{b} = 44880.87 ~(+3.66)$$
 
"""

# ╔═╡ a1dca343-7819-4ed2-9639-067805c17085
md"""
# Results: $$\nu_R\nu_L\beta\beta$$
ROI: $$\varphi \in (5, 175)^{\circ}, E_{sum} \in (100, 3000)keV, r \in (0, 50)mm$$. **Basically take everything?**

 $$T^{1/2} \geq 1.30 \times 10^{22}$$ y; $$\varepsilon = 0.08$$; $$\bar{b} = 115520.50$$
"""

# ╔═╡ 8fc9a87c-ec29-4028-a20f-9c8cdec4c866
md"""
# Conclusions

- Using n-dimensional approach (sometimes) leads to *slight* increase in sensitivities
  - there are issues that need investigation
  - better optimization should improve the results even further
     - The choice of data-cut on $$r<50 mm$$ **is way too strict!**
     - Will try to optimize ToF values as well
  - Choice of variables can be extended, at the cost of computing time
  - Machine learning techniques for MVA could improve it even further
- Bayesian approach is more sensitive than frequentist when signal shape is more different from backgrounds!
  - not sure what's wrong with RH spectra, but this could be it?
- Need to implement neutron data
- We have an idea to use data-partitioning for angular distribution


"""

# ╔═╡ 9e7dbfe3-a38b-45ef-b6c6-74ca9c3392df
md"""
# Future
Based on research of angular distribution for the electrons, depending on which side of foil they are emitted, the reconstruction is differently efficient. We could use this fact!
![same_side](https://github.com/Shoram444/SNSensitivityEstimate/blob/main/notebooks/Sensitivity_exotic_Bayes_nDim/same_side_opposite.png?raw=true)
"""

# ╔═╡ 8789bf21-3712-404c-b047-b0304d7917da


# ╔═╡ 1f4ae42b-61be-44ed-bb82-60e7cfb582ba
begin
	files_directory = "fal5_8perc_Boff_TKrec_evis_bcu_J38"

	data_processes = load_data_processes(
	    files_directory, 
	    "sumE",
	    fwhm = 0.0
	)
	signal = get_process("bb0nu_foil_bulk", data_processes) |> first

    background = [
        get_process("bb_foil_bulk", data_processes) |> first,
        get_process("Bi214_foil_bulk", data_processes) |> first,
        get_process("Bi214_wire_surface", data_processes) |> first,
        get_process("Tl208_foil_bulk", data_processes) |> first,
        get_process("K40_foil_bulk", data_processes) |> first,
        get_process("Pa234m_foil_bulk", data_processes) |> first,
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
    sig_hist = restrict(get_bkg_counts_1D(signal), Emin, Emax)


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
        lines!(ax, midpoints(binedges(sig_hist)), (bincounts(sig_hist) .* (1e24/Thalf)) .+ bincounts(sum(bkg_hists)), label = L"0\nu\beta\beta", color = :red, linewidth = 3.5)

        ax.yticks = ([1e-5, 1e-3, 1e-1, 1e1, 1e3, 1e5], [L"10^{-5}",L"10^{-3}", L"10^{-1}", L"10^{1}", L"10^{3}", L"10^{5}"])
        ax.xticks = 0:500:3500
        Legend(f[2,1], ax, orientation=:horizontal, fontsize=8, nbanks = 3)
        f
    end
end

# ╔═╡ da4b1d6e-3447-4451-b5b8-fbff4ebfe3ee
let
	f = Figure()
	a1 = Axis(f[1,1], xlabel = "energy", ylabel = "a.u.", title = "signal vs background")
	hsig = Hist1D(signal.dataVector; binedges = signal.bins) |> normalize
	h1 = stairs!(a1, hsig, label = "signal")



	hbkg = Hist1D(background[1].dataVector; binedges = background[1].bins) |> normalize
	hbkg.bincounts .*= 10
	h2 = stairs!(a1, hbkg, label = "background")

	xband_min = 2700
	xband_max = 3100
	ymin = 0
	ymax = 1.2maximum(bincounts(hbkg))
	poly!(
		a1,
		[xband_min, xband_max, xband_max, xband_min],
    	[ymin, ymin, ymax, ymax],
		color = (:blue, 0.3),
		label = "ROI"
	)
	axislegend(a1, position = :lt)
	
	f
	
end

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

# ╔═╡ c8d16dc1-b06a-4b61-b7cc-0304dadfbcb8
begin
	vars = [
	    "phi", 
	    "sumE"
	]
	bins = (
	    phi = (0,180),
	    sumE = (0, 3500)
	)
	processesND = load_ndim_processes("fal5_TKrec", bins, vars)
	nothing
end

# ╔═╡ 4cb301ac-87ba-46a9-969d-8bc2ae77bdc4
begin
	signalND = get_process("bb0nu_foil_bulk", processesND) |> first
	backgroundND = get_process("Tl208_foil_bulk", processesND) |> first
	backgroundNDbb = get_process("bb_foil_bulk", processesND) |> first
	set_signal!(backgroundNDbb, false)
	set_nTotalSim!( signalND, 1e8 )
	set_nTotalSim!( backgroundND, 1e8 )
	set_nTotalSim!( backgroundNDbb, 1e8 )
	nothing
end

# ╔═╡ e842c11a-81db-47f1-aa65-4d53ee4bfea5
begin
	phiSignalND = getproperty.(signalND.data,:phi)
	esumSignalND = getproperty.(signalND.data,:sumE)
	phiBackgroundND = vcat(getproperty.(backgroundND.data,:phi), getproperty.(backgroundNDbb.data,:phi))
	esumBackgroundND = vcat(getproperty.(backgroundND.data,:sumE), getproperty.(backgroundNDbb.data,:sumE))
	nothing
end

# ╔═╡ fef5daa1-3806-450d-a80f-e3416f066924
let 
	f = Figure(size=(1300,600))
	a1 = Axis(f[1,1], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Signal")
	
	h2dsignal = Hist2D(
			(phiSignalND, esumSignalND);
			binedges = (0:5:180, 0:100:3500)
		)
	h2dbackground = Hist2D(
			(phiBackgroundND, esumBackgroundND);
			binedges = (0:5:180, 0:100:3500)
		)

	h1 = plot!(a1, h2dsignal, colormap = :coolwarm)
	text!(p_l, E_l-300, text="ROI", fontsize = 15, color = :red)
	

	a2 = Axis(f[1,2], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Background")
	h2 = plot!(a2, h2dbackground, colormap = :coolwarm)
	text!(p_l, E_l-300, text="ROI", fontsize = 15, color = :red)


	poly!(
		a1,
		[p_l,p_u,p_u,p_l],
		[E_l,E_l,E_u,E_u],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	poly!(
		a2,
		[p_l,p_u,p_u,p_l],
		[E_l,E_l,E_u,E_u],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	f
end

# ╔═╡ f913fc9a-ebec-4686-9cbb-c8cbac70c26b
let 
	signalNDM1 = get_process("bb0nuM1_foil_bulk", processesND) |> first
	set_nTotalSim!( signalNDM1, 1e8 )
	phiSignalNDM1 = getproperty.(signalNDM1.data,:phi)
	esumSignalNDM1 = getproperty.(signalNDM1.data,:sumE)

	
	f = Figure(size=(1300,800))
	a01 = Axis(f[1,1:2], xlabel = L"E_{sum}", title = L"$E_{sum}$ $0\nu\beta\beta\chi^0$")
	a02 = Axis(f[1,3:4], xlabel = L"\varphi", title = L"$E_{sum}$ $0\nu\beta\beta\chi^0$")
	h1signalE = Hist1D(esumSignalNDM1; binedges = (0:100:3500))
	h1signalPhi = Hist1D(phiSignalNDM1; binedges = (0:5:180))

	plot!(a01, h1signalE)
	plot!(a02, h1signalPhi)
	
	a1 = Axis(f[2,1], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Signal efficiency")
	
	h2dsignal = Hist2D(
			(phiSignalNDM1, esumSignalNDM1);
			binedges = (0:5:180, 0:100:3500)
		) 
	h2dsignal.bincounts .*= inv(1e8)
	h2dbackground = Hist2D(
			(phiBackgroundND, esumBackgroundND);
			binedges = (0:5:180, 0:100:3500)
		) |> normalize
	h2dbackground.bincounts .*= 1e5

	h1 = plot!(a1, h2dsignal, colormap = :coolwarm)
	Colorbar(f[2,2], h1)

	a2 = Axis(f[2,3], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Background counts")
	h2 = plot!(a2, h2dbackground, colormap = :viridis)
	Colorbar(f[2,4], h2)


	poly!(
		a1,
		[10,175,175,10],
		[2500,2500,3100,3100],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	poly!(
		a2,
		[10,175,175,10],
		[2500,2500,3100,3100],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	
	
	f
end

# ╔═╡ 5583f1a0-8383-42d8-88fb-e3d6e806954c
let 
	signalNDM2 = get_process("bb0nuM2_foil_bulk", processesND) |> first
	set_nTotalSim!( signalNDM2, 1e8 )
	phiSignalNDM2 = getproperty.(signalNDM2.data,:phi)
	esumSignalNDM2 = getproperty.(signalNDM2.data,:sumE)

	
	f = Figure(size=(1300,800))
	a01 = Axis(f[1,1:2], xlabel = L"E_{sum}", title = L"$E_{sum}$ $0\nu\beta\beta\chi^0\chi^0$")
	a02 = Axis(f[1,3:4], xlabel = L"\varphi", title = L"$E_{sum}$ $0\nu\beta\beta\chi^0\chi^0$")
	h1signalE = Hist1D(esumSignalNDM2; binedges = (0:100:3500))
	h1signalPhi = Hist1D(phiSignalNDM2; binedges = (0:5:180))

	plot!(a01, h1signalE)
	plot!(a02, h1signalPhi)
	
	a1 = Axis(f[2,1], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Signal efficiency")
	
	h2dsignal = Hist2D(
			(phiSignalNDM2, esumSignalNDM2);
			binedges = (0:5:180, 0:100:3500)
		) 
	h2dsignal.bincounts .*= inv(1e8)
	h2dbackground = Hist2D(
			(phiBackgroundND, esumBackgroundND);
			binedges = (0:5:180, 0:100:3500)
		) |> normalize
	h2dbackground.bincounts .*= 1e5

	h1 = plot!(a1, h2dsignal, colormap = :coolwarm)
	Colorbar(f[2,2], h1)

	a2 = Axis(f[2,3], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Background counts")
	h2 = plot!(a2, h2dbackground, colormap = :viridis)
	Colorbar(f[2,4], h2)


	poly!(
		a1,
		[0,180,180,0],
		[1200,1200,3100,3100],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	poly!(
		a2,
		[0,180,180,0],
		[1200,1200,3100,3100],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	
	
	f
end

# ╔═╡ 87f06108-c73d-4509-8e97-4f1ea07dd484
let 
	signalNDdRH = get_process("RH037_foil_bulk", processesND) |> first
	set_nTotalSim!( signalNDdRH, 1e8 )
	phiSignalNDdRH = getproperty.(signalNDdRH.data,:phi)
	esumSignalNDdRH = getproperty.(signalNDdRH.data,:sumE)

	
	f = Figure(size=(1300,800))
	a01 = Axis(f[1,1:2], xlabel = L"E_{sum}", title = L"$E_{sum}$ $\nu_R\nu_L\beta\beta$")
	a02 = Axis(f[1,3:4], xlabel = L"\varphi", title = L"$E_{sum}$ $\nu_R\nu_L\beta\beta$")
	h1signalE = Hist1D(esumSignalNDdRH; binedges = (0:100:3500))
	h1signalPhi = Hist1D(phiSignalNDdRH; binedges = (0:5:180))

	plot!(a01, h1signalE)
	plot!(a02, h1signalPhi)
	
	a1 = Axis(f[2,1], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Signal efficiency")
	
	h2dsignal = Hist2D(
			(phiSignalNDdRH, esumSignalNDdRH);
			binedges = (0:5:180, 0:100:3500)
		) 
	h2dsignal.bincounts .*= inv(1e8)
	h2dbackground = Hist2D(
			(phiBackgroundND, esumBackgroundND);
			binedges = (0:5:180, 0:100:3500)
		) |> normalize
	h2dbackground.bincounts .*= 1e5

	h1 = plot!(a1, h2dsignal, colormap = :coolwarm)
	Colorbar(f[2,2], h1)

	a2 = Axis(f[2,3], xlabel = L"\varphi", ylabel = L"E_{sum}", title = "Background counts")
	h2 = plot!(a2, h2dbackground, colormap = :viridis)
	Colorbar(f[2,4], h2)


	poly!(
		a1,
		[0,180,180,0],
		[100,100,3000,3000],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	poly!(
		a2,
		[0,180,180,0],
		[100,100,3000,3000],
		color = (:black,0),
		strokewidth = 3,
		strokecolor = :red,
	 )
	
	
	f
end

# ╔═╡ Cell order:
# ╟─1a041346-47c7-41a4-a30c-5710a5dd8786
# ╟─33db2852-a205-4a2c-8744-037f5b7a6a80
# ╟─107241e3-a8b7-43af-b26f-b2b205f7e294
# ╟─838099af-cd31-480e-8952-9eda9df9e10a
# ╟─65b97922-2625-43a4-bf9f-622cce900bcf
# ╟─4011c3d8-d298-40f1-9986-77e991a25152
# ╟─781740e1-5f0c-4b99-834c-1ac6a75118b7
# ╟─4db07780-248f-4658-82a5-507ea65edb0b
# ╟─16fbb9bb-d86b-433b-91fa-a12f8008eaa9
# ╟─5099a914-3f46-4e86-9570-bc13f5d3ca80
# ╟─de64766d-fda7-4a0e-bc7c-80bb90554bac
# ╟─bad82399-f4f0-401c-acc7-b7f67155f27f
# ╟─c12114cf-0b28-4579-8bbb-147e8acc6f33
# ╟─4a2ca349-25ff-4a16-a4cf-075767700240
# ╟─eb22cfcd-a08e-4e57-9f02-113612804fc1
# ╟─aae87c04-2352-4db5-9018-f8d0fa75e4b9
# ╟─c8539419-4919-48a9-ad41-0dcba22318e8
# ╟─6edb8cab-331c-403c-8c93-b8af84f38989
# ╟─da4b1d6e-3447-4451-b5b8-fbff4ebfe3ee
# ╟─1e127dad-951d-4873-b711-fb8a23f3d712
# ╟─c519cbff-afb1-4ac7-9bfe-1e07cf78979b
# ╟─755fcdeb-d7a8-4d41-bb7d-2ec9b5e3065b
# ╟─d7354883-f530-4bfe-955d-ca37263d0c7b
# ╟─87c50992-3347-47e8-8cb4-3f0dcee41ebf
# ╟─9ae6e7cf-fa54-4d17-834f-a11f7201bdff
# ╟─112a6be9-e6ac-4ee2-b7b4-82fc7d5d7b88
# ╟─fef5daa1-3806-450d-a80f-e3416f066924
# ╟─8bc26abf-c593-44d8-b1d1-5f2db8d07705
# ╠═4cb301ac-87ba-46a9-969d-8bc2ae77bdc4
# ╠═e842c11a-81db-47f1-aa65-4d53ee4bfea5
# ╟─6fa76137-73fc-47ea-927b-0ffed51149f7
# ╟─cdcc36a8-888b-42e4-8721-18e183a3ae5e
# ╟─724e4ff0-6dea-4a73-8708-0c886c98a8db
# ╟─8680f3a5-6ed9-40de-b8a3-c2e34b4507e0
# ╟─e6bbd64e-a2eb-462f-837a-079ef07f3c43
# ╟─e768fb37-96d1-4181-9e36-1aa5af23db1f
# ╟─8c041449-ec1c-4bdf-ac7b-0ef59ba777a4
# ╟─1910efa8-1530-4139-8287-d9d05c3e69ae
# ╟─5b254b0e-6b4b-46be-bed9-8df0efa3de44
# ╟─e83967c4-4bf9-4715-8aca-d022209c77d5
# ╟─518693c6-694b-45eb-a08b-c3a449a0d0bf
# ╟─ea51b339-1d49-4f04-a8fe-744f9e855026
# ╟─ed6b7268-412d-45a9-895c-7248b0686efe
# ╟─c57ec8c2-68d1-47b9-b7e7-0a415ce85031
# ╟─9043a542-5102-488e-a3ab-0f377648145a
# ╟─e5664e79-6d50-43dc-963a-d2207c48b9a8
# ╟─cc0d49b8-3d6b-4a2c-8300-0aea834df028
# ╟─25f1adc4-a6f4-4a10-9fe1-b0d1817f0ece
# ╟─085a4b0a-1fdc-44d8-a8ce-7158ac54256f
# ╟─187dfcf3-48d8-4ee9-8d33-8b3bac063b37
# ╟─8a0763e0-bf30-49df-930c-d28ded622ded
# ╟─6d3bd21c-a856-407e-b3e7-d98d8cad0c9b
# ╟─559c3769-c956-47e7-b5ea-07b866333104
# ╟─4e8c3fa4-b9b4-4656-a61e-2756ce10f3f1
# ╟─bad14d75-c550-4221-884c-d57c787f3039
# ╟─b2a28a93-2524-43e3-a9fa-6ef30ac213af
# ╟─8e5fa765-c8a9-4419-9e4f-c58a213be563
# ╟─917eca81-62a3-49e5-a9d1-47d2b6f9c129
# ╟─a1d02f00-354a-4ebb-b3f0-388c2a69c210
# ╟─7da00dac-d0a8-4b77-973e-9e46e5fa964d
# ╟─eb626008-caed-4c21-b563-1f0db822355f
# ╟─f913fc9a-ebec-4686-9cbb-c8cbac70c26b
# ╟─e2e4a637-dc3e-4cd5-a156-9cd0ebe56495
# ╟─5583f1a0-8383-42d8-88fb-e3d6e806954c
# ╟─a1dca343-7819-4ed2-9639-067805c17085
# ╟─87f06108-c73d-4509-8e97-4f1ea07dd484
# ╟─8fc9a87c-ec29-4028-a20f-9c8cdec4c866
# ╟─9e7dbfe3-a38b-45ef-b6c6-74ca9c3392df
# ╟─96f0e67c-16af-11f0-08a5-ffefdd2c8429
# ╟─39a6343a-b95c-4ea4-98e2-83eab6d6c602
# ╟─c38eeb3e-11c2-40d0-98cb-89711b799441
# ╟─992a6366-86e0-4d00-ab9c-ac29190ff693
# ╟─7d6e5825-c187-4d55-b877-1e1193190ff2
# ╟─8789bf21-3712-404c-b047-b0304d7917da
# ╟─1f4ae42b-61be-44ed-bb82-60e7cfb582ba
# ╟─c8d16dc1-b06a-4b61-b7cc-0304dadfbcb8
