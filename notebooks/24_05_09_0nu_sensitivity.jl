### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e0d5731a-086d-11ef-0315-5705b559453d
using DrWatson

# ╔═╡ 2c1253c6-5601-4701-944c-03fcad58c175
# ╠═╡ show_logs = false
@quickactivate "SNSensitivityEstimate"

# ╔═╡ 773e7c01-dc0f-464a-9696-9a89e2302aa8
# ╠═╡ show_logs = false
begin
	using CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Measurements, PrettyTables, Turing, ReverseDiff
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
	include(scriptsdir("LoadData.jl"))
end

# ╔═╡ a5b84247-8464-4a80-af07-5309d7bcaef9
using PlutoUI

# ╔═╡ e4c0246e-b67a-4776-8f15-6efd5aec1f00
html"<button onclick='present()'>present</button>"

# ╔═╡ 3bd341f6-eae9-41f5-b5eb-b6c3e76cbaf2
html"""
<style>
  main {
    max-width: 1500px;
	font-size: 22px;
  }
</style>
"""

# ╔═╡ c56b1bac-c440-4f71-a7a3-088b3042dfad
dao_theme = dao();

# ╔═╡ 099b4434-ca9d-4748-ba3e-7f097e3229e3
set_theme!(dao_theme);

# ╔═╡ 4b9a7ce2-3e5f-4995-af22-985b5fb0f14e
md"""
# Sensitivity to $0\nu\beta\beta$ in SuperNEMO: two ways! Frequentist vs Bayessian approach. 

### Also compared to Cupid-0
### Maros Petro
### SuperNEMO Analysis meeting 9.5.2024

"""

# ╔═╡ acdec082-7a61-489f-9642-9a459ab9f04a
md"""
# Outline

- **Introduction to sensitivity**
- **Simulated data-set**
- **Frequentist approach: "The Feldman Cousins approach"**
- **Bayessian approach: "Here's what the data tells me..."**
"""

# ╔═╡ bee405ae-e2bd-48b8-8de8-652b99bbd8c4
md"""
# 1. Introduction to sensitivity
"""

# ╔═╡ 955810a1-9ef6-4531-8156-02b070c29676
md"""
# But what is sensitivity really? 

In a "simple" experiment trying to estimate the half-life of a rare decay we can use the following forumla:

```math 
T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{n_{OBS}} 
```

!!! danger "Issues with this equation:"
	- There is always background present (even in background-less experiments)
	- Radioactive decay is a statistical process - even very low probability is **not** zero probability
	- What if the number of observed **signal counts** is $n_{OBS} = 0$? Infinite half-life? wtf..
"""

# ╔═╡ 11b39aa3-2c3b-4d5f-8ef0-e938c3cb3399
md"""
### Statistics to the rescue! (It is both our friend and our enemy...)

- In case of no observation, we turn the problem on its head: 
> Instead of asking **what is the half-life of the process which we haven't measured** we ask the following question: **What is the lowest value of half-life we are able to measure with our detector?**

**Mathematically speaking we do the following adjustments to the half-life formula:**

1. We exchange $n_{OBS}$ with an **interval of values** which represents our *belief* that $n_{OBS}$ lies within the limits of: 
```math
\mu_L \leq n_{OBS} \leq \mu_U
```

2. Since $T^{1/2}$ is inversly proportional to $n_{OBS}$ we get:

```math
ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_U} \leq T^{1/2} \leq ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_L} ~~~ \rightarrow ~~~ T^{1/2}\in(T^{1/2}_{lower}, T^{1/2}_{upper})
```

3. Since we expect
```math 
n_{OBS} = 0 \rightarrow \mu_L = 0 \rightarrow T^{1/2}_{up} = \infty
```

4. We are left with **only the lower bound** for the half-life $\Rightarrow$ *sensitivity* 
```math
T^{1/2} \geq ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_U} 
```
"""

# ╔═╡ e4c6c2b0-2afe-4667-894c-71c8476239c6
md"""
# Now about the $\mu_U$... 

!!! info "There's a few issues here:"
	1. In a *real* experiment, there is background. How do we incorporate it to the calulation? 
	2. **Everything is statistical!** (Fluctuations of background and signal measurements...) 
	3. What about the chosen confidence interval (CI)? Central, one-sided, other? The choice of CI **affects** calculated $\mu_U$. But that is dangerous grounds!

#### We need 3 ingredients:
1. Background model
2. Choice of Confidence level - let's go with 90% (standard)
3. Methodology for calculating $\mu_U$

### Which approach do we chose?
#### Frequentist vs Bayess? ...Both!

"""

# ╔═╡ a65e5feb-fe4a-41a0-9cc9-530c0c90a317
md"""
# 2. Simulated data-set
"""

# ╔═╡ 21b41886-2a18-40d2-b57c-62862b30b8fe
md"""
# To demonstrate how the sensitivity is calculated, I show a small example of calculating the SuperNEMO's sensitivity to $0\nu\beta\beta$ 

#### Obviously, first we must have data. We can simulate the **signal** and **relevant background processes** to get the estimate of what we can expect. Here's the simulation set-up:


1. Falaise version 5.1.1 
2. Reconstruction pipeline: SimRC, MockCalibration, **CAT**, TrackFit, ChargedParticleTracker
3. Magnetic field on (default 25G)
4. Standard energy resolution of 8% at 1MeV (for Main wall OMs)
5. Foil geometry: `RealisticFlat`
6. Signal: `Se82.0nubb`
7. Background: `Se82.2nubb`, `Bi214_foil_bulk`, `Bi214_foil_surface`, `Bi214_wire_surface`, `Tl208_foil_bulk`


"""

# ╔═╡ cede3473-fae8-4998-ba6e-01f9d05ad54c
md"""
## Here's what the raw data looks like:
"""

# ╔═╡ 6821c511-ef71-401e-9f36-1719886e7b54
begin
	signalProcessName = "bb0nu_foil_bulk"
	sigProcessESum = get_process(signalProcessName, probProcessesESum)
	
	bb2nuProcess = get_process("bb_foil_bulk", probProcessesESum)
	bb2nuProcess.signal = false
	bkgProcessesESum = [
	    bb2nuProcess,
	    get_process("Bi214_foil_bulk", probProcessesESum),
	    get_process("Bi214_foil_surface", probProcessesESum),
	    get_process("Bi214_wire_surface", probProcessesESum),
	    get_process("Tl208_foil_bulk", probProcessesESum),
	]
	nothing
end

# ╔═╡ 52165d49-4a46-41d5-9803-20047b1ce12e
h1d_signal_raw = Hist1D( sigProcessESum.dataVector; binedges= binningDict[:SumE]);

# ╔═╡ 3cbcfaf5-b367-40e2-9093-d10f852d403a
h1d_background_raw = [ Hist1D(p.dataVector;binedges= binningDict[:SumE]) for p in bkgProcessesESum  ];

# ╔═╡ f57a6ed8-2332-4ab1-abef-27c44881c7af
md"""
## First we take a look at the stacked and *normalized* histograms of all processes together. 

#####  This will help in seeing the shapes of the spectra and their "reach". 
"""

# ╔═╡ 84c47066-c7b7-4653-8e7a-9d24305fe091
# ╠═╡ show_logs = false
let
	f = Figure(size = (1000, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "stacked normalized histogram")

	labels= ["signal", "2nubb", "Bi214_foil_bulk", "Bi214_foil_surface", "Bi214_wire_surface", "Tl208_foil_bulk"]

	st = FHist.stackedhist!(ax, vcat(normalize(h1d_signal_raw), normalize.(h1d_background_raw)), color= Makie.wong_colors())
	elements = [PolyElement(polycolor = st.attributes.color[][i]) for i in 1:length(labels)]
	f[1,2] = axislegend(ax, elements, labels,)
	f
end

# ╔═╡ ec440471-bd5a-44ca-bd8c-7be723c07179
md"""
!!! note "About the figure:"
	- The peak of 0nu (in blue) is smeared quite a bit due to energy resolution + electron energy losses in gas. *In an ideal detector this should be delta function.*
	- Furthermore, we see that each of the processes depicted pose as background for 0nu.
"""

# ╔═╡ c06c24af-9387-4094-b8f5-ac89e01230c0
md"""
# Now let's take a look at the estimated spectrum we should see after 2.88yr of measurement.

### To estimate the spectrum we use the following equation:
```math
n = A \cdot m \cdot \varepsilon \cdot t
```

- Basically to calculate the number of expected events for each process we need to know the activity, time of measurement, mass or volume, efficiency. 
!!! update "What simulation can provide!"
	The only value we get from simulation is the **efficiency**! All the other variables in the equation are experimental and cannot be simulated!
"""

# ╔═╡ ce3b0244-4acf-4670-bfca-9d9715ee0d93
md"""
> All the parameters I use in this analysis are listed [here](https://github.com/Shoram444/SNSensitivityEstimate/blob/main/scripts/Params.jl)
"""

# ╔═╡ 8d9a138e-e015-411c-bc27-966d06fdb901
md"""
## Estimated background only spectrum for SuperNEMO experiment: 17.5kgyr. 

!!! note "Disclaimer!"
	- I'm only showing background processes relevant to 0nu. There are others present at lower energies. 
"""

# ╔═╡ 54824dab-3172-46f2-b130-31f55a8262fb
h1d_background = get_bkg_counts_1D.(bkgProcessesESum);

# ╔═╡ 3c5fadf4-5791-4f33-9ce1-6394611b977b
# ╠═╡ show_logs = false
let
	f = Figure(size = (1000, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1:2,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "stacked backgrounds histogram",xticklabelrotation=45)
	ax2 = Axis(f[2,2], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "zoomed to 2500,3500", limits=(2500,3500, nothing,nothing), yaxisposition = :right, xticklabelrotation=45)
	
	labels= ["2nubb", "Bi214_foil_bulk", "Bi214_foil_surface", "Bi214_wire_surface", "Tl208_foil_bulk"]
	st = FHist.stackedhist!(ax, h1d_background, color= Makie.wong_colors())
	st = FHist.stackedhist!(ax2, h1d_background, color= Makie.wong_colors())

	ylims!(ax2, 0, 5)
	elements = [PolyElement(polycolor = st.attributes.color[][i]) for i in 1:length(labels)]
	f[1,2] = axislegend(ax, elements, labels,)
	colgap!(f.layout, 1, Relative(0.07))
	f
end

# ╔═╡ bceec9c1-785e-4aec-8195-82d1db4ee838
# used to set ylims
min_cts = minimum(@. minimum( filter(x-> x>0, bincounts(h1d_background)) ) );

# ╔═╡ f5f2503a-1809-49c2-8839-3ea17c99b286
md"""
# We can take a look at the same plot but in log scale:
"""

# ╔═╡ 84acc905-d854-42c2-a862-e362410205b2
let
	f = Figure(size = (1000, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1:2,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "stacked backgrounds histogram \nLog-scale", yscale =log10, xticklabelrotation=45)
	ax2 = Axis(f[2,2], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "zoomed in", yscale =log10, limits=(2500,3500, nothing,nothing), yaxisposition = :right, xticklabelrotation=45)
	
	labels= ["2nubb", "Bi214_foil_bulk", "Bi214_foil_surface", "Bi214_wire_surface", "Tl208_foil_bulk"]
	st = hist!(ax, sum(h1d_background), label =labels[1],c=1)
	for i=2:length(h1d_background)
		hist!(ax, sum(h1d_background[i:end]), label=labels[i], c=i)
	end

	st2 = hist!(ax2, sum(h1d_background), label =labels[1],c=1)
	for i=2:length(h1d_background)
		hist!(ax2, sum(h1d_background[i:end]), label=labels[i], c=i)
	end
	
	ylims!(ax, min_cts, 1e5)
	ylims!(ax2, min_cts, 1e1)
	f[1, 2] = Legend(f, ax)
	colgap!(f.layout, 1, Relative(0.07))
	f
end

# ╔═╡ 78c184cf-e6c2-485d-9876-e3811ae00c8f
md"""
# 3. Frequentist approach: "The Feldman Cousins approach"
"""

# ╔═╡ 6ed9cdf8-6d18-456e-934f-035f8e2df1d3
md"""
# Why do we care about background anyway? 

### Because in a frequentist approach $\mu_U$ depends **only on background**: $\mu_U \Rightarrow \mathcal S(\bar{b})$

- where $\bar{b}$ is *expected number of background counts in ROI* (in units counts).
"""

# ╔═╡ 883a548d-54ec-4242-9e65-38118045a517
md"""
#### We introduce a statistical factor $\mathcal S(\bar{b})$ which represents *limiting value on background counts, above which the experiment can claim discovery*.

```math
\mathcal S(\bar{b}) \begin{cases} 
      F.C. & \mathrm{if ~\bar{b}~ is~ small} \\
      Gauss & \mathrm{if ~\bar{b}~ is~ large} \\
   \end{cases}
```

Detailed explanation for what $\mathcal S(\bar{b})$ means is given on [1](https://theses.hal.science/tel-02122290), [2](https://theses.cz/id/334te1/STAG97246.pdf). You may know this under a different name, $N_{excl}$: *number of excluded bkg events*. 

For now let's just say it's a number we calculate and replace $\mu_U$ with in the sensitivity equation. 



"""

# ╔═╡ 9be29518-9b98-4f6d-a504-40a27ada1107
md"""
# Skipping forward a few steps, which I already presented on in [5833-v3](http://nile.hep.utexas.edu/cgi-bin/DocDB/ut-nemo/private/ShowDocument?docid=5833)...

#### For each combination of ROI we can calculate the sensitivity map. The best sensitivity is obtained as a maximum of the map:
"""

# ╔═╡ bb31751d-eaa8-4b7e-aad0-cf894a6ee7f7
α=1.64;

# ╔═╡ 7e9ac44f-b1af-477b-a9ab-762e12a9d294
# ╠═╡ show_logs = false
t12MapESum = get_tHalf_map(SNparams, α, sigProcessESum, bkgProcessesESum...);

# ╔═╡ e5cba772-d3a1-4b5b-9897-164183366e04
let
	f = Figure( size = (800, 700), fontsize= 22, )
	ax = Axis(f[1,1], xlabel="min_ROI [keV]", ylabel="max_ROI [keV]", title = "sensitivity map as a function of ROI", aspect=1)
	p = plot!(ax, t12MapESum, colorscale=log10)
	Colorbar(f[1,2], p, label="log sensitivity [yr]", scale=log10)
	colsize!(f.layout, 1, Aspect(1, 1.0))
	f
end

# ╔═╡ 196c8517-cb7a-4286-b630-98488d8bbe18
best_t12ESum = get_max_bin(t12MapESum);

# ╔═╡ abc26356-b76e-4bc0-a1fd-c9aed18c0841
expBkgESum = get_bkg_counts_ROI(best_t12ESum, bkgProcessesESum...);

# ╔═╡ be30be54-ccad-4e7e-a758-f55b4dfbfa35
effbb = lookup(sigProcessESum, best_t12ESum);

# ╔═╡ a4b6a249-d490-4c97-af08-aaa7548d0b86
nexcl = get_FC(effbb, α);

# ╔═╡ 6879f6a2-6d14-45de-99d6-f0bf0baceb85
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3);

# ╔═╡ a859ecc0-bfae-4c5b-87da-5d5804a60ba1
b_in_ROI = [get_bkg_counts_ROI(best_t12ESum, p) for p in bkgProcessesESum];

# ╔═╡ 010d6d22-df4c-4dbf-b31e-21629201b8c6
bkg_isotope_names = [p.isotopeName for p in bkgProcessesESum];

# ╔═╡ cd9e7bcf-193c-4c2d-9e42-5fd9f0004886
b_eff_in_ROI = [ lookup(p, best_t12ESum)*100 for p in bkgProcessesESum];

# ╔═╡ ff1cc95d-bb3f-49b0-b988-8cdeca763b9d
df_iso_bkg_ROI = DataFrame(isotope = bkg_isotope_names, b = b_in_ROI, efficiency = b_eff_in_ROI);

# ╔═╡ 6356676f-32a7-452d-b6ee-3d6267c9e9fa
pretty_isotopoes=pretty_table(HTML,df_iso_bkg_ROI, header=["isotope", "b", "eff [%]"]);

# ╔═╡ 14d69080-2f1b-48a6-8c05-d58b5b4d3992
md"""
## From this map we can find the following: 

- best ROI: ($( best_t12ESum[:minBinEdge] ) - $( best_t12ESum[:maxBinEdge] )) keV 
- signal efficiency in ROI: $(round(effbb, sigdigits=3) *100) %
- expected background counts $\bar{b}$ in ROI: $(round(expBkgESum, sigdigits=3) ) 
- corresponding $\mathcal S(\bar{b})$: $(round(nexcl, sigdigits=3) ) 
- ### max sensitivity: $(round(best_t12ESum[:maxBinCount], sigdigits=3)) yr

### Background overview table:
$(
pretty_isotopoes
)
"""

# ╔═╡ 46086903-ec96-4f86-8081-ddfb713fde30
md"""
# 4. Bayessian approach: "Here's what the data tells me..."

!!! danger "Danger zone!" 
	❗The following section is based on my (so far) **limited**  understaing of Bayessian analysis. ❗

	There's still much to learn for me, so please take the results with a grain of 🧂.
"""

# ╔═╡ 05c42c02-7306-412b-b205-74701de21634
md"""
## Bayes Theorem
"""

# ╔═╡ e54eb759-9274-4df3-9cc3-12ba8925b067
md"""
## Frequentist vs Bayess statistics:

- Frequency --> Belief
- Point-estimate --> distribution
- Parameters fixed --> parameters distributed over pdf
- Confidence interval --> credibility interval
"""

# ╔═╡ 339077ec-a94d-4d36-9003-3dc690cdc00a
md"""
## MCMC sampling
"""

# ╔═╡ d941baf4-d60a-40f8-a35d-bf76e16dbeda
md"""
## Bayess approach to $\mu_U$: 

- In frequentist approach we exchanged $\mu_U$ for a background dependent value calculated from FC approach
- In bayessian approach we calculate $\mu_U$ **directly**!! (Or at least our belief in what it should be...)
"""

# ╔═╡ 742f7fd5-01e7-4731-a7a1-ad84bb86b8c7
md"""
## Bayess approach to SuperNEMO Sensitivity. An example:

!!! note "Methodology"
	We start again with calculation of half-life:

	$T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{n_{OBS}}$

	Now remember, we made the change where: $\mu_L \leq n_{OBS} \leq \mu_U$
	- In Frequentist approach we made a **point-estimate** on the value of $\mu_U$ by calculating $\bar{b}$. (We have $\mu_U = \mathcal S(\bar{b})$.)
	- In Bayessian approach, we consider $n_{OBS}$ as a parameter of the fit and **infer** the underlying distribution by MCMC fit. 
	- The output of the MCMC fit is the **distribution** of $n_{OBS}$, from which we can directly calculate the credible interval: $(\mu_L, \mu_U)$!!! 🎆
	

"""

# ╔═╡ 278988c0-6da6-46f2-a50d-89ea32edb69a
md"""
## Applying model to data...
"""

# ╔═╡ 4cd71a16-3bb8-4d81-b1f7-f67eb6cad196
md"""

!!! update "Likelihood, prior, posterior...:"
	As desribed above. The Bayessian analysis takes advantage of Bayess Theorem. This means we must define **likelihood and prior** to calculate **posterior**. 
	Therefore we must make a few assumptions in making our **model**: 
	1. What is the likelihood (what do we believe is the underlying distribution of our measured data)?
	2. What are our fit parameters and what are their priors?

	To build the model we look at the data:

"""

# ╔═╡ 28a5b0b3-8055-49a4-8623-688c785be07d

let
f = Figure(size = (600, 400), fontsize= 16, fonts = (; regular = "TeX"))
ax = Axis(f[1:2,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "stacked backgrounds histogram", limits=(2500,3500, nothing,nothing))

labels= ["2nubb", "Bi214_foil_bulk", "Bi214_foil_surface", "Bi214_wire_surface", "Tl208_foil_bulk"]
st = FHist.stackedhist!(ax, h1d_background, color= Makie.wong_colors())

ylims!(ax, 0, 1.1*lookup(sum(h1d_background), 2500))
elements = [PolyElement(polycolor = st.attributes.color[][i]) for i in 1:length(labels)]
f[1,2] = axislegend(ax, elements, labels,)
f
end


# ╔═╡ 699d1d0b-e9c0-4ac7-bfc0-05b36fd3ab23
md"""
> I chose a ROI (2500-3500 keV) slightly larger than the ROI predicted by Frequentist approach (2700-3100 keV). 
"""

# ╔═╡ 2c04e398-8c31-4f9b-a00a-528c9bf60358
md"""

!!! update "Likelihood, prior, posterior...:"
	As with any attempts at making a model with the hopes of describing data, we must make a number of assumption:
	1. Data consists of some amount of background $n_b$ and some amount of signal $n_s$ (i.e. the Signal + Background model)
	2. Background is distributed exponentially with decay constant $\lambda$: $pdf(b) = \lambda e^{-\lambda E}$ 

"""

# ╔═╡ 273f495a-02b8-494f-80c8-77d7ec83a6ab


# ╔═╡ 0dc66ff5-c2c5-4b86-a971-f9f01e7dc7ff
md"""
# The end! 😏
"""

# ╔═╡ aa75ac59-2e73-440e-90ed-ff397c696892
@bind mu PlutoUI.Slider(5:15)

# ╔═╡ 5a7f967c-32f4-448e-8382-39547104a499
f, ax, p = plot(Normal(mu),);

# ╔═╡ 62a82632-ef0d-4617-8b3d-981ab675d56a
xlims!(ax, 0, 100)

# ╔═╡ 7d7e0703-12d1-4f37-8998-f966a7af6e62
f

# ╔═╡ Cell order:
# ╟─e4c0246e-b67a-4776-8f15-6efd5aec1f00
# ╟─3bd341f6-eae9-41f5-b5eb-b6c3e76cbaf2
# ╟─e0d5731a-086d-11ef-0315-5705b559453d
# ╟─2c1253c6-5601-4701-944c-03fcad58c175
# ╟─c56b1bac-c440-4f71-a7a3-088b3042dfad
# ╟─099b4434-ca9d-4748-ba3e-7f097e3229e3
# ╟─773e7c01-dc0f-464a-9696-9a89e2302aa8
# ╟─4b9a7ce2-3e5f-4995-af22-985b5fb0f14e
# ╟─acdec082-7a61-489f-9642-9a459ab9f04a
# ╟─bee405ae-e2bd-48b8-8de8-652b99bbd8c4
# ╠═955810a1-9ef6-4531-8156-02b070c29676
# ╟─11b39aa3-2c3b-4d5f-8ef0-e938c3cb3399
# ╟─e4c6c2b0-2afe-4667-894c-71c8476239c6
# ╟─a65e5feb-fe4a-41a0-9cc9-530c0c90a317
# ╟─21b41886-2a18-40d2-b57c-62862b30b8fe
# ╟─cede3473-fae8-4998-ba6e-01f9d05ad54c
# ╟─6821c511-ef71-401e-9f36-1719886e7b54
# ╟─52165d49-4a46-41d5-9803-20047b1ce12e
# ╟─3cbcfaf5-b367-40e2-9093-d10f852d403a
# ╟─f57a6ed8-2332-4ab1-abef-27c44881c7af
# ╟─84c47066-c7b7-4653-8e7a-9d24305fe091
# ╟─ec440471-bd5a-44ca-bd8c-7be723c07179
# ╟─c06c24af-9387-4094-b8f5-ac89e01230c0
# ╟─ce3b0244-4acf-4670-bfca-9d9715ee0d93
# ╟─8d9a138e-e015-411c-bc27-966d06fdb901
# ╟─3c5fadf4-5791-4f33-9ce1-6394611b977b
# ╟─54824dab-3172-46f2-b130-31f55a8262fb
# ╟─bceec9c1-785e-4aec-8195-82d1db4ee838
# ╟─f5f2503a-1809-49c2-8839-3ea17c99b286
# ╟─84acc905-d854-42c2-a862-e362410205b2
# ╟─78c184cf-e6c2-485d-9876-e3811ae00c8f
# ╟─6ed9cdf8-6d18-456e-934f-035f8e2df1d3
# ╟─883a548d-54ec-4242-9e65-38118045a517
# ╟─9be29518-9b98-4f6d-a504-40a27ada1107
# ╟─e5cba772-d3a1-4b5b-9897-164183366e04
# ╟─7e9ac44f-b1af-477b-a9ab-762e12a9d294
# ╟─bb31751d-eaa8-4b7e-aad0-cf894a6ee7f7
# ╟─14d69080-2f1b-48a6-8c05-d58b5b4d3992
# ╟─196c8517-cb7a-4286-b630-98488d8bbe18
# ╟─abc26356-b76e-4bc0-a1fd-c9aed18c0841
# ╟─be30be54-ccad-4e7e-a758-f55b4dfbfa35
# ╟─a4b6a249-d490-4c97-af08-aaa7548d0b86
# ╟─6879f6a2-6d14-45de-99d6-f0bf0baceb85
# ╟─a859ecc0-bfae-4c5b-87da-5d5804a60ba1
# ╟─010d6d22-df4c-4dbf-b31e-21629201b8c6
# ╟─cd9e7bcf-193c-4c2d-9e42-5fd9f0004886
# ╟─ff1cc95d-bb3f-49b0-b988-8cdeca763b9d
# ╟─6356676f-32a7-452d-b6ee-3d6267c9e9fa
# ╟─46086903-ec96-4f86-8081-ddfb713fde30
# ╠═05c42c02-7306-412b-b205-74701de21634
# ╠═e54eb759-9274-4df3-9cc3-12ba8925b067
# ╠═339077ec-a94d-4d36-9003-3dc690cdc00a
# ╠═d941baf4-d60a-40f8-a35d-bf76e16dbeda
# ╟─742f7fd5-01e7-4731-a7a1-ad84bb86b8c7
# ╟─278988c0-6da6-46f2-a50d-89ea32edb69a
# ╠═4cd71a16-3bb8-4d81-b1f7-f67eb6cad196
# ╟─28a5b0b3-8055-49a4-8623-688c785be07d
# ╟─699d1d0b-e9c0-4ac7-bfc0-05b36fd3ab23
# ╠═2c04e398-8c31-4f9b-a00a-528c9bf60358
# ╠═273f495a-02b8-494f-80c8-77d7ec83a6ab
# ╠═0dc66ff5-c2c5-4b86-a971-f9f01e7dc7ff
# ╠═a5b84247-8464-4a80-af07-5309d7bcaef9
# ╠═5a7f967c-32f4-448e-8382-39547104a499
# ╠═62a82632-ef0d-4617-8b3d-981ab675d56a
# ╠═aa75ac59-2e73-440e-90ed-ff397c696892
# ╠═7d7e0703-12d1-4f37-8998-f966a7af6e62
