### A Pluto.jl notebook ###
# v0.19.41

#> [frontmatter]
#> title = "Sensitivity to 0nubb in SuperNEMO: two ways! Frequentist vs Bayessian approach."
#> 
#>     [[frontmatter.author]]
#>     name = "Maros Petro"

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
	using CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Measurements, PrettyTables, Turing, ReverseDiff, LaTeXStrings
	push!(LOAD_PATH, srcdir())
	using SensitivityModule
	include(scriptsdir("Params.jl"))
	nothing
end

# ╔═╡ d33cdccb-3c3e-4394-b16b-6a0bd56b5a3f
using BAT, DensityInterface, IntervalSets, BinnedModels

# ╔═╡ a5b84247-8464-4a80-af07-5309d7bcaef9
using PlutoUI

# ╔═╡ 752369c9-464f-44fc-97df-0f642bc83af7
using Optim

# ╔═╡ e4c0246e-b67a-4776-8f15-6efd5aec1f00
html"<button onclick='present()'>present</button>"

# ╔═╡ 3bd341f6-eae9-41f5-b5eb-b6c3e76cbaf2
html"""
<style>
  main {
    max-width: 1500px;
	font-size: 150%;
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
### SuperNEMO Analysis meeting 16.5.2024

"""

# ╔═╡ acdec082-7a61-489f-9642-9a459ab9f04a
md"""
# Outline

#### 1. A general comparison of Bayessian and frequentist statistics
#### 2. Introduction to sensitivity
#### 3. Simulated data-set 
#### 4. Model
#### 5. Frequentist approach: "The Feldman Cousins approach"
#### 6. Bayessian approach: "Here's what the data tells me..."
"""

# ╔═╡ 503673d7-7dac-43eb-a493-0aea57674720
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 50e7af06-2bce-4914-a867-ec6261a6714b
md"""
#### There's a few things that distinguish the two approaches to statistics: 

!!! update " "Philosophycal" differences:" 
	**How we ask about hypotheses...**
	- *Frequentist*: "Under given hypothesis, how unlikely is the observed result?" ➡ **p-value**
	- *Bayessian*: "Given the observed results, how likely is the hypothesis?" ➡ **posterior distribution**

	**How we treat probability...**
	- *Frequentist*: Probability is long-run frequency!
	- *Bayessian*: Probability is a degree of belief!
"""

# ╔═╡ 839d78ac-2ad9-4cc0-92df-7d3402a5b291
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 29274571-9856-4dd7-9afa-d5204fac46ca
md"""
#### There's a few things that distinguish the two approaches to statistics: 

!!! note " "Practical" differences:" 
	**How we treat parameters...**
	- *Frequentist*: True value of a parameter is fixed and unknown
	- *Bayessian*: Parameters are treated as random variables with their own distributions 

	**How we estimate parameters...**
	- *Frequentist*: Point estimates (i.e. MLE); uncertainty must be estimated
	- *Bayessian*: Posterior distribution (updated belief of prior); uncertainty included in distribution

	**How we treat intervals...**
	- *Frequentist*: **confidence interval** constructed so that in an ensemble of experimenents $\alpha$ fraction contains parameter
	- *Bayessian*: **credible interval** represents the range within which the parameter lies with a certain probability ($\alpha$), given the observed data and prior information
"""

# ╔═╡ 50b657c7-f326-4cfa-a4dc-160980d1ec24
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 80b557d9-96b8-4d00-93cd-2986798a6054
md"""
#### Bayess Theorem: 

```math
P( H | E ) = \frac{P(E|H)\cdot P(H)}{P(E)}
```

where:
-  $P(H∣E)$: Conditional/Posterior probability of hypothesis H given evidence E
-  $P(E∣H)$: Likelihood of evidence E given hypothesis H
-  $P(H)$: Prior probability of hypothesis H
-  $P(E)$: Marginal likelihood of evidence E

#### Bayess theorem describes the probability of an event, based on prior knowledge of conditions that might be related to the event.
"""

# ╔═╡ 93d6b82b-6d0e-4be5-812b-5cee189987f2
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 7a925206-47b3-4702-9df2-4074239ca4e6
md"""
#### Bayess Theorem: Not just for Bayessians statistics...

Bayess theorem is not exclusive to bayessian statistics. An (amazing) example use of application of Bayess theorem in frequentist framework for COVID tests by [3Blue1Brown](https://www.youtube.com/watch?v=lG4VkPoG3ko)! 

#### The devil's in the details. (Specifically how Bayess theorem is used in Bayessian statistics.)

"""

# ╔═╡ ce031aaf-ceb5-4b82-8e0e-bc3a795daecd
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 50493bc0-eb22-4af5-ac47-5da71eb96632
md"""

In Bayessian statistics (inference) one attemtps to **infer** the parameter's posterior distribution given **data** and **priors**. **We update the knowledge of prior probability!**
- From Bayess theorem: $posterior \propto likelihood \times prior$ 
- In some cases we can ignore marginal probability (it's generally difficult to calculate but really serves as just a normalization constant)
- In bayessian inference we do the following "not-so-formal" change to Bayess Theorem:
```math
	p(\theta | data) \propto p(data | \theta) \times p(\theta)
```
Now:
-  $p(\theta | data)$ is the **updated** posterior probability **distribution** of parameter $\theta$ **given** data
-  $p( data | \theta )$ is the likelihood **distribution** (a distribution that maps all plausible ways data can be obtained)
-  $p(\theta)$ is the prior **distribution** of parameters - generally a **subjective** prior as this is up to the choice of analyzer. 

"""

# ╔═╡ 74d7f5c7-d5bb-41c1-9e72-b5da45cb3538
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ cbf182c3-4e2e-445d-9814-f5f6146a25e4
md"""
### Bayessian process:
1. Start with a prior distribution which reflects the **subjective** belief of how the parameters can be distributed
2. Define likelihood of data given parameter
3. Collect data and apply Bayess theorem to obtain updated posterior
4. **Posterior becomes the new prior**
5. Repeat with new data...

"""

# ╔═╡ 7f8276ae-0a48-40ad-bdec-ee2cefc7ec55
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 7be26605-60c2-4a03-a090-385caa299091
md"""
!!! danger "Sampling the distribution"
	In general, it is very difficult to sample the full distributions. Using a random sampling processes like inverse CDF or rejection methods are computationally very costly!
"""

# ╔═╡ 68da87d9-c5bb-4aa5-837b-6803fb7687e6
md"""
### Markov Chain Monte Carlo (MCMC) in Bayesian Inference:

The Bayessian process described in previous slide is generally implemented via the use of the MCMC methods:
> In statistics, Markov chain Monte Carlo (MCMC) is a class of algorithms used to draw samples from a probability distribution. Given a probability distribution, one can construct a Markov chain whose elements' distribution approximates it – that is, the Markov chain's equilibrium distribution matches the target distribution. The more steps that are included, the more closely the distribution of the sample matches the actual desired distribution. (source: Wikipedia)

The general idea of a markov chain is that the **next element drawn learns from the previously drawn element**. 
"""

# ╔═╡ ee33bb4c-3f6b-4be3-be77-18e121208444
md"""
# 1. A general comparison of Bayessian and frequentist statistics
"""

# ╔═╡ 3a34adbb-bc3c-4043-a524-e0697523900d
md"""
### Summary:

- Bayessian statistics: 
  * **A statistical process that updates the probability for a hypothesis as more evidence or information becomes available.**
- Bayess theorem:
  * It can be used in frequentist approach.
  * In Bayessian statistics, use probability distributions - continuosly update.
- Bayessian inference process: 
  1. Define model: likelihood and prior,
  2. Gather data,
  3. Update posterior,
  4. Repeat. 
- MCMC: 
  * A class of algorithms that **efficienctly** sample from distributions. 
"""

# ╔═╡ bee405ae-e2bd-48b8-8de8-652b99bbd8c4
md"""
# 2. Introduction to sensitivity
"""

# ╔═╡ 955810a1-9ef6-4531-8156-02b070c29676
md"""
# But what is sensitivity really? 

Imagine an ideal experiment with no background whatsoever. In this "simple" experiment we are trying to estimate the half-life of a rare decay. To do so we can use the following forumla:

```math 
T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{n_{S}} 
```

!!! danger "Issues with this equation:"
	- If the process is allowed (remember, $0\nu\beta\beta$ could really be forbidden) we expect to measure some number of events $n$. 
	-  $n$ does not necessarily have to equal the **expected number of signal events** $n_S$.
	- Radioactive decay is a statistical process - we may be lucky and get many events, or we may be unlucky and get nothin'...
	- What do we do if the number of observed **counts** is $n = 0$? Infinite half-life? 🤔
"""

# ╔═╡ 11b39aa3-2c3b-4d5f-8ef0-e938c3cb3399
md"""
# Statistics to the rescue! (It is both our friend and our enemy...)

!!! update "In case of no signal observed, we turn the problem on its head: "
	Instead of asking **what is the half-life of the process which we haven't measured** we ask the following question: **What is the lowest value of half-life we are able to measure with our detector?**

**Mathematically speaking we do the following adjustments to the half-life formula:**

1. We exchange $n_{S}$ with an **interval of values** which represents our *belief* that $n_{S}$ lies within the limits of: 
```math
\mu_L \leq n_{S} \leq \mu_U
```

2. Since $T^{1/2}$ is inversly proportional to $n_{S}$ we get:

```math
ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_U} \leq T^{1/2} \leq ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_L} ~~~ \rightarrow ~~~ T^{1/2}\in(T^{1/2}_{lower}, T^{1/2}_{upper})
```

3. Since we expect
```math 
n_{S} = 0 \rightarrow \mu_L = 0 \rightarrow T^{1/2}_{up} = \infty
```

4. We are left with **only the lower bound** for the half-life $\Rightarrow$ *sensitivity* 
```math
T^{1/2} \geq ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{\mu_U} 
```
"""

# ╔═╡ e4c6c2b0-2afe-4667-894c-71c8476239c6
md"""
# Now about the $\mu_U$ (and $\mu_L$)... 

#### There are two ways of looking at what it means:
1. ## Bayessian approach:
!!! note "Bayessian terminology:" 
	- In Bayessian analysis, $n_S$ is a parameter of intereset and is distributed according to some pdf $p(n_S | data)$. 
	- A **credible interbal** is a range of plausible values that $n_S$ can take, i.e. a 90% credible interval covers 90% of plausible values. (The choice of central, one-sided, highest density is up to the scientist. **And it matters!**) 
	
	 **$(\mu_L, \mu_U)$ are then the lower and upper boundaries of the 90% credible interval:** $\int_{\mu_L}^{\mu_U} p(n_S|data) dn_S = 90 \%$.

"""

# ╔═╡ 94ad3e49-34d8-45e5-ac65-a9241136f226
let 
	f = Figure(fontsize = 24, fontfamily="TeX")
	ax = Axis(f[1,1], xlabel= L"n_S", title = "90% credible interval on \nPosterior distribution of signal counts")
	clr = "#56B4E9"
	p = plot!(ax, Erlang(3,0.4), color = clr)

	xs = range(0.0, 2.7, 100)
	y_l =zeros(100)
	y_u = pdf.(Erlang(3,0.4), xs)
	band!(ax,xs, y_l, y_u, color =(clr, 0.3))
	hideydecorations!(ax, minorgrid=:false, )
	hidexdecorations!(ax, label = false, ticklabels = false)
	hidespines!(ax, :l, :t, :r)
	ax.xticks= ([0, 2.7], [L"\mu_L", L"\mu_U"])
	xlims!(ax, 0,7)
	f
end

# ╔═╡ 875130d7-f0e3-4be0-a963-c427712ccdcf
md"""
2. ## Frequentist:
!!! note "Frequentist terminology:"
	- In frequentist approach, the parameter $n_S$ is unknown and **fixed**. 
	- A **confidence interval** is built in such a way that if an ensemble of identical experiments were repeated many times, a fraction - say 90% - of the calculated intervals would cover $n_S$. 
	 **$\mu_U$ is the upper bound of the confidence interval - it doesn't actually claim probability of $n_S$ being within it!** (We would have to repeat experiment many times for that!)
"""

# ╔═╡ 557e5713-6563-4e4e-b0ce-a9ee67cf12e6
let
function generate_data(n::Int, true_mean::Real, true_std::Real)
    data = rand(Normal(true_mean, true_std), n)
    return data
end

function calculate_confidence_intervals(data::Vector{Float64}, confidence_level::Float64)
    n = length(data)
    mean_val = mean(data)
    std_err = std(data) / sqrt(n)
    z_critical = quantile(Normal(), (1 + confidence_level) / 2)

    lower_bound = mean_val - z_critical * std_err
    upper_bound = mean_val + z_critical * std_err

    return lower_bound, upper_bound
end

function plot_confidence_intervals(true_mean::Real, true_std::Real, num_samples::Int, confidence_level::Float64)
    
    intervals = []
    for _ in 1:num_samples
		data = generate_data(num_samples, true_mean, true_std)
        lower, upper = calculate_confidence_intervals(data, confidence_level)
        push!(intervals, (lower, upper))
    end

	x = 1:num_samples
    y = [(intervals[i][1] + intervals[i][2]) / 2 for i in 1:num_samples]
    errors = [(intervals[i][2] - intervals[i][1]) / 2 for i in 1:num_samples]

    fig = Figure(fontsize = 24, fontfamily="TeX")
    ax = Axis(fig[1, 1])

    p = scatter!(ax, y, x, color=:black, markersize=3)
    for i in 1:num_samples
		if( intervals[i][1] > true_mean || intervals[i][2] < true_mean )
        	linesegments!(ax, [(intervals[i][1], x[i]), (intervals[i][2], x[i])], color=:red, linewidth=4)
		else
        	linesegments!(ax, [(intervals[i][1], x[i]), (intervals[i][2], x[i])], color="#56B4E9", linewidth=4)
		end
    end
    errorbars!(ax, y, x, errors, color=:blue)

    vlines!(ax, [true_mean], color=:red, linestyle=:dash, label="True Mean")
    ax.title = "Confidence Intervals (90% cover true n_S)"
    ax.ylabel = "Sample Number"
    ax.xlabel = L"n_S"
	hideydecorations!(ax, minorgrid=:false, )
	hidexdecorations!(ax, label = false, ticklabels = false)
	hidespines!(ax, :l, :t, :r)
	ax.xticks= ([5], ["true"])
	xlims!(ax, 4.4, 5.6)
	
    fig, ax, p
end

f,a, p = plot_confidence_intervals(5, 1, 40, 0.90)
f
end

# ╔═╡ b1991653-1363-4810-a78d-923c4ce64daa
md"""

# Furthermore real experiments have background! If we measure $n$ counts, we may expect that some of that originates from *background*, $n_B$, and (hopefully) some from *signal*, $n_S$. 

!!! info "There's a few issues here...again:"
	1. How do we incorporate background into the calculation?
	2. How do we calculate confidence/credible intervals with background? 

> Today, I will show an example of calculating sensitivity to $0\nu\beta\beta$ for SuperNEMO using both approaches, Bayessian and frequentist. I will do so on a generated *sample data-set*. 

#### However, today I will only show how this is done in **Bayessian analysis**. The details of the frequentist method will be presented at some other time in the future. (If you want, you can read on it in Docdb: [4816-v1](http://nile.hep.utexas.edu/cgi-bin/DocDB/ut-nemo/private/ShowDocument?docid=4816).)

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
7. Background 500M events each: `Se82.2nubb`, `Bi214_foil_bulk`, `Bi214_foil_surface`, `Bi214_wire_surface`, `Tl208_foil_bulk`


"""

# ╔═╡ cede3473-fae8-4998-ba6e-01f9d05ad54c
md"""
# Here's what the raw data looks like:
"""

# ╔═╡ 6821c511-ef71-401e-9f36-1719886e7b54
# ╠═╡ show_logs = false
begin
	allProcesses = load_processes("fal5_8perc", "sumE")
	signalProcessName = "bb0nu_foil_bulk"
	sigProcessESum = get_process(signalProcessName, allProcesses)
	
	bb2nuProcess = get_process("bb_foil_bulk", allProcesses)
	bb2nuProcess.signal = false
	bkgProcessesESum = [
	    bb2nuProcess,
	    get_process("Bi214_foil_bulk", allProcesses),
	    get_process("Bi214_foil_surface", allProcesses),
	    get_process("Bi214_wire_surface", allProcesses),
	    get_process("Tl208_foil_bulk", allProcesses),
	]
	nothing
end

# ╔═╡ 52165d49-4a46-41d5-9803-20047b1ce12e
h1d_signal_raw = Hist1D( sigProcessESum.dataVector ; binedges= binningDict[:SumE]);

# ╔═╡ 3cbcfaf5-b367-40e2-9093-d10f852d403a
h1d_background_raw = [ Hist1D(p.dataVector  ;binedges= binningDict[:SumE] ) for p in bkgProcessesESum  ];

# ╔═╡ f57a6ed8-2332-4ab1-abef-27c44881c7af
md"""
## First we take a look at the stacked and *normalized* histograms of all processes together. 

#####  This will help in seeing the shapes of the spectra and their "reach". 
"""

# ╔═╡ 84c47066-c7b7-4653-8e7a-9d24305fe091
# ╠═╡ show_logs = false
let
	f = Figure(size = (1000, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1:2,1], xlabel = "energy [keV]", ylabel = "normalized counts / 100keV" , title = "stacked normalized histogram", xticklabelrotation=45)
	ax2 = Axis(f[2,2], xlabel = "energy [keV]", ylabel = "" , title = "zoomed to 2500,3500", limits=(2500,3500, nothing,nothing), yaxisposition = :right, xticklabelrotation=45)

	labels= ["signal", "2nubb", "Bi214_foil_bulk", "Bi214_foil_surface", "Bi214_wire_surface", "Tl208_foil_bulk"]

	st = FHist.stackedhist!(ax, vcat(normalize(h1d_signal_raw), normalize.(h1d_background_raw)), color= Makie.wong_colors())
	st2 = FHist.stackedhist!(ax2, vcat(normalize(h1d_signal_raw), normalize.(h1d_background_raw)), color= Makie.wong_colors())
	elements = [PolyElement(polycolor = st.attributes.color[][i]) for i in 1:length(labels)]
	f[1,2] = axislegend(ax, elements, labels,)
	colgap!(f.layout, 1, Relative(0.07))
	f
end

# ╔═╡ ec440471-bd5a-44ca-bd8c-7be723c07179
md"""
!!! note "About the figure:"
	- The peak of $0\nu\beta\beta$ (in blue) is smeared quite a bit due to energy resolution + electron energy losses in gas. *In an ideal detector this should be delta function.*
	- Furthermore, we see that each of the processes depicted pose as background for 0nu.
"""

# ╔═╡ c06c24af-9387-4094-b8f5-ac89e01230c0
md"""
# Now let's take a look at the estimated spectrum we should see after 2.88yr of measurement. (This corresponds to 17.5kgy exposure.)

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
# Estimated background only spectrum for SuperNEMO experiment: 17.5kgyr. 


> Disclaimer: I'm only showing background processes relevant to 0nu. There are others present at lower energies. 
"""

# ╔═╡ 54824dab-3172-46f2-b130-31f55a8262fb
h1d_background = get_bkg_counts_1D.(bkgProcessesESum);

# ╔═╡ 3c5fadf4-5791-4f33-9ce1-6394611b977b
# ╠═╡ show_logs = false
let
	f = Figure(size = (1000, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1:2,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "stacked backgrounds histogram \nlinear scale",xticklabelrotation=45)
	ax2 = Axis(f[2,2], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "zoomed to 2500,3500", limits=(2500,3500, nothing,nothing), yaxisposition = :right, xticklabelrotation=45)
	
	labels= ["2nubb", "Bi214_foil_bulk", "Bi214_foil_surface", "Bi214_wire_surface", "Tl208_foil_bulk"]
	st = FHist.stackedhist!(ax, h1d_background, color= Makie.wong_colors())
	st = FHist.stackedhist!(ax2, h1d_background, color= Makie.wong_colors())

	ylims!(ax2, 0, 10)
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

# ╔═╡ 68231ef5-9425-4e4c-98a8-586323938d45
md"""

"""

# ╔═╡ 1c44a17c-41a5-4ec2-a37e-b4c20e9407cc


# ╔═╡ 0033c567-5cf4-4067-bc50-00543dabc798


# ╔═╡ c69bf1d5-e7ae-4b31-8fe3-cbb872e68424
md"""
# To make the analysis fair when compring the two approaches I use the same sample data:

!!! update "Gererating sample spectrum from Background model"
	The spectra above show the bacgkround model - *the relative background rates per bin*. From this model we can generate a sample spectrum in the following way:
	- Take the bin height in the above spectrum - this represents the expected background rate in 17.5kgy exposure per keV, $\bar{b_i}$.
	- For each bin $i$, generate a random sample of $n_i$ background events Poisson distributed with $\mu_i = \bar{b_i}$. This is can be done with `rand(Poisson(mu))`
	- Dataset then consists of $n_i$ uniformaly distributed energies within each $b_i$
	- Combine to form a sample spectrum.

I am generating background only events, as I assume $n_S = 0$. This will be explained in more detail later.
"""


# ╔═╡ 58241414-5ec2-463a-8422-f7bd8e30b606
md"""
## Generating data...
"""

# ╔═╡ 91546e37-2d22-4616-9d19-a4d973e80e54
@bind clicked PlutoUI.Button("Generate data")

# ╔═╡ d5b67056-4212-40b6-b25c-007cacd0509a
md"""
!!! note "Note"
	Since this type of data generation is a statistical process, it fluctuates a lot! 
"""

# ╔═╡ 06e1024a-5afb-4cd8-a75e-4bcfb4a5ab8f
h1d_full_bkg = sum(h1d_background);

# ╔═╡ a47c9336-33de-4c89-bc64-3e1e12036243
function generate_sample_data(h1d::Hist1D)
	n_i = map(x->rand(Poisson(x)), bincounts(h1d))
	b_i = collect(binedges(h1d))
	bin_widths = diff(b_i) # vector of bin widths (in case of non-uniform binning)
	
	data = Vector{Float64}(undef, sum(n_i))
	idx_slice = 1
	for (ni, bi, bw_i) in zip(n_i, b_i, bin_widths)
		if(ni != 0)
			data[idx_slice:idx_slice+ni-1] = rand(Uniform(bi, bi + bw_i ), ni)
			idx_slice+=ni
		end
	end
	return data 
end;

# ╔═╡ 873c693a-23a9-4ec2-ac3a-04cd9c9ac40b
md"""
# From now on we will use this sample set to perform both analyses on. 
"""

# ╔═╡ 52e4668e-9804-492a-9deb-830f8cdedfcf
md"""
# 3. Model
"""

# ╔═╡ fab2761b-fa89-40e1-bfd8-e30d5524b7c2
md"""
!!! note "Analysis procedure"
	In both the frequentist and bayessian approach, our goal is to estimate the parameters that describe our measured data. The general idea is:
	1. Get data
	2. Create model
	3. Fit model to data
	4. Extract parameters of interest.
"""

# ╔═╡ f768a9e3-7ab7-4c36-bf5b-2737a638c325
md"""
# 3. Model

### Now that we have data, we must create the model. 

For the frequentist approach, we will use Maximum Likelihood Estimate to fit the data, which means we must define the **likelihood**: 

#### There are two choices we can make here. A detailed model and a general model.
   1. **The detailed model** would be to build the likelihood from the simulated spectra of **each individual background (and signal) component** and define a *mixture model*:
```math
	\mathcal L( x | \vec{\theta} ) = \frac{1}{\sum \theta_i}\prod_i^{N^{obs}}\left( \theta_1 * pdf(sig, E_i) + \sum_{i=2} \theta_i * pdf(bkg_i, E_i) \right)
```
Here, $N^{obs}$ are observed events. $pdf(sig, E_i)$ and $pdf(bkg_i, E_i)$ are the signal and background normalized spectra evaluated at energies $E_i$. The parameters of interest $\vec\theta$ are the activities. (In case of background only model $\theta_1 = 0$.) 
	
  

"""

# ╔═╡ 370440df-155a-46c2-a3f1-75f4ebf73bc9
md"""
# 3. Model
"""

# ╔═╡ 255a4f62-964b-47d6-96f0-663cfcdcf09f
md"""
 2. **The general model** is to build the likelihood as the sum of **gaussian** distributed signal and **exponentially** distributed **total** background:
```math
	\mathcal L( x | n_S, n_B, \lambda ) = \frac{1}{n_S+n_B}\prod_i^{N^{obs}}{\left( n_B\lambda e^{-\lambda E_i} + n_S \frac{1}{\sigma\sqrt{2\pi}}e^{-(\frac{E_i-Q}{\sigma})^2}\right)}
```
Here, the parameters of interest are:
-  $n_S$ expected signal counts (in case of no-background model $n_S = 0$)
-  $n_B$ expected total background counts
-  $\lambda$ total background decay rate

#### I will use the general model for both bayessian and frequentist analyses.
"""

# ╔═╡ 78c184cf-e6c2-485d-9876-e3811ae00c8f
md"""
# 4. Frequentist approach: "Brought to you by Feldman and Cousins"
"""

# ╔═╡ 6ed9cdf8-6d18-456e-934f-035f8e2df1d3
md"""
# We come back to the question, why do we care about background anyway? 

### Because in a frequentist approach $\mu_U$ depends **only on background**: $\mu_U \Rightarrow \mathcal S(\bar{b})$

- where $\bar{b}$ is *expected number of background counts in ROI* (in units counts).
"""

# ╔═╡ 883a548d-54ec-4242-9e65-38118045a517
md"""
#### We introduce a statistical factor $\mathcal S(\bar{b})$ which represents *limiting value on background counts, above which the experiment can claim discovery*.

```math
\mathcal S(\bar{b}) = \begin{cases} 
      F.C. & \mathrm{if ~\bar{b}~ is~ small} \\
      Gauss & \mathrm{if ~\bar{b}~ is~ large} \\
   \end{cases}
```

!!! note "Details on FC"
	Detailed explanation for what $\mathcal S(\bar{b})$ means is given on [1](https://theses.hal.science/tel-02122290), [2](https://theses.cz/id/334te1/STAG97246.pdf). 
	You may know this under a different name, $N_{excl}$: *number of excluded bkg events*. 

### For now let's just say it's a number we calculate and replace $\mu_U$ with in the sensitivity equation. 



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
	Colorbar(f[1,2], p, label="sensitivity [yr]", scale=log10)
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
nexcl = get_FC(expBkgESum, α);

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
# From this map we can find the following values for the background model: 

- ### best ROI: ($( best_t12ESum[:minBinEdge] ) - $( best_t12ESum[:maxBinEdge] )) keV 
- signal efficiency in ROI: $(round(effbb, sigdigits=4) *100) %
- expected background counts $\bar{b}$ in ROI: $(round(expBkgESum, sigdigits=3) ) 
- corresponding $\mathcal S(\bar{b})$: $(round(nexcl, sigdigits=3) ) 
- ### best sensitivity: $(round(best_t12ESum[:maxBinCount], sigdigits=3)) yr

### ROI Background overview table:
$(
pretty_isotopoes
)
> Radon is by far the largest bkg contributor...and I'm using a very optimistic activity.

!!! danger "High expected background?"
	In comparison to the "standard" ``\bar{b}`` = 0.5 cts and the one I calculated earlier this year ``\bar{b} \approx 1.0`` cts, the new value ``\bar{b} = 2.55`` cts is **way too big**. Two things have changes that may be the culprit:
	1. I used `SimRC` conditions (with realistic **flat** foil though!).
	2. I generated a **lot more events** which decreased relative error in high energy regions.

"""

# ╔═╡ 6ea1a763-3268-4f77-83a0-f6065d77eb0a
md"""
# However, in order to be able to **compare** the sensitivity with Bayessian analysis, we have to calculate the sensitivity **NOT** from the background model but from the **sample data-set**. We do this in the following way: 

!!! note "Frequentist procedure for data-set"
	- First we cut the sample data-set to just energies close to the ROI (region a bit wider that ROI): $E \in (a,b)$.
	- This data-set we use as a fitting region. I use the **general model** defined above with MLE fit. 
	- We integrate the model over the boundaries of the best ROI to obtain $\bar{b}$.
	- We convert sample's $\bar{b}$ to $\mathcal S(\bar{b})$.
	- Lastly, plug into sensitivity equation and calulate sample's sensitvitiy.
"""

# ╔═╡ 69f34564-f7e5-45d5-90d6-501591ae67ae
md"""
# Generate sample data, choice of fitting region, fit...
"""

# ╔═╡ 2485d73a-28bc-461a-be00-0083dba8a5fb
@bind clicked1 PlutoUI.Button("Generate data")

# ╔═╡ 796d7a86-dae3-4e0c-9dc7-157d245425e9
md"""
a = $(@bind a PlutoUI.Slider(1600:100:2700, default=2300; show_value=true)) keV

b = $(@bind b PlutoUI.Slider(3100:100:3500, default=3400; show_value=:true)) keV
"""

# ╔═╡ 46086903-ec96-4f86-8081-ddfb713fde30
md"""
# 5. Bayessian approach: "Here's what the data tells me..."

!!! danger "Danger zone!" 
	❗The following section is based on my (so far) **limited**  understaing of Bayessian analysis. ❗

	There's still much to learn for me, so please take the results with a grain of 🧂.
"""

# ╔═╡ d941baf4-d60a-40f8-a35d-bf76e16dbeda
md"""
### Bayess approach to $\mu_U$: 

- In frequentist approach we exchanged $\mu_U$ for a background dependent value calculated from FC approach
- In bayessian approach we calculate $\mu_U$ **directly**!! (Or at least our belief in what it should be...)
"""

# ╔═╡ 742f7fd5-01e7-4731-a7a1-ad84bb86b8c7
md"""
# Bayess approach to SuperNEMO Sensitivity. An example:

!!! note "Methodology"
	We start again with calculation of half-life:

	$T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{n_{S}}$

	Now remember, we made the change where: $\mu_L \leq n_{S} \leq \mu_U$
	- In Frequentist approach we made a **point-estimate** on the value of $\mu_U$ by calculating $\bar{b}$. (We have $\mu_U = \mathcal S(\bar{b})$.)
	- In Bayessian approach, we consider $n_{S}$ as a parameter of the fit and **infer** the underlying distribution by applying MCMC sampling and employing Bayess theorem. 
	- The output of the procedure is the **distribution** of $n_{S}$, from which we can directly calculate the credible interval: $(\mu_L, \mu_U)$!!! 🎆
	

"""

# ╔═╡ 278988c0-6da6-46f2-a50d-89ea32edb69a
md"""
# Applying model to data...
"""

# ╔═╡ 4cd71a16-3bb8-4d81-b1f7-f67eb6cad196
md"""

!!! update "Likelihood, prior, posterior...:"
	As desribed above. The Bayessian analysis takes advantage of Bayess Theorem. This means we must define **likelihood and prior** to calculate **posterior**. 
	Therefore we must make a few assumptions in making our **model**: 
	1. What is the likelihood (what do we believe is the underlying distribution of our measured data)? Here, we will again use the **general model** likelihood defined previously. 
	2. What are our **prior beliefs** we put on the parameters?
"""

# ╔═╡ 14ef85c1-4be6-45c3-b2f8-cae5d5a21aa7
md"""
# Defining the model:
"""

# ╔═╡ 2c04e398-8c31-4f9b-a00a-528c9bf60358
md"""

!!! update "Likelihood, prior, posterior...:"
	We make a *choice* to use the **general model** likelihood:
	```math
		\mathcal L(data | n_S, n_B, \lambda) = \frac{1}{n_S+n_B}\prod_i^{N^{obs}}{\left( n_B\lambda e^{-\lambda E_i} + n_S \frac{1}{\sigma\sqrt{2\pi}}e^{-(\frac{E_i-\mu}{\sigma})^2}\right)}
	```

	The parameters of interest, now **❗distributions❗** we want to *infer* from the data are: 
	-  $p(\lambda|data)$ the background decay constant
	-  $p(n_B|data)$ the number of expected background events
	-  $p(n_S|data)$ the number of expected signal events

	The (somewhat informative) priors are:
	-  $p(\lambda|data) \sim Uniform(10^{-3}, 10^3)$
	-  $p(n_B|data) \sim Uniform( 0.0, 1000.0 )$ 
	-  $p(n_S|data) \sim Uniform( 0.0, 10.0 )$*

	*- the prior for $n_S$ is bounded by the "possible values" calculated from $n = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{T^{1/2}}$. Where for the lower and upper bound, I subsituted $T^{1/2} = \infty$ and $T^{1/2} = 10^{24} y$ (order of magnitude of the best limit at the moment)

"""

# ╔═╡ 6376b61a-f6a2-4dbd-9ba1-b2bcd96f5d9d
md"""
# Sampling from posterior:
"""

# ╔═╡ 2134a77e-34bb-4b7f-a268-23fd1e6d52a1
md"""
I will sample the posterior within the MCMC chain using the No-U-Turn (NUTS) algorithm with 65% acceptance rate, with 2 chains and 5000 samples per chain. 
"""

# ╔═╡ b9047089-49d5-4ad7-8e2c-d2029176b393
@bind train PlutoUI.Button("Sample Now!!") 

# ╔═╡ 7ca9d29b-688b-45a3-bb13-dce16e067fae
md"""
!!! update " "
	The output of the sampling process is a `chains` object which holds the information on the individual samples and their statistics:
"""

# ╔═╡ 14bde55d-7695-431c-930f-2b7181665937
md"""
# Visualizing the chains to validate the output:
"""

# ╔═╡ dcc61c0b-5529-42b7-8fd1-88033b060b0b
md"""
# Visualizing the model fit to data to validate the output:
"""

# ╔═╡ 9a8cff73-af96-4c5b-a455-4b890f9c3dc3
md"""
#### The plot below shows the fit using Bayessian inference along with various CIs'.  
"""

# ╔═╡ c2bdf038-475e-4b52-8c4a-92e1ad09c1e0
md"""
!!! danger " "
	One of the main advantages here is the straight-forward evaluation of fit uncertainty! 
"""

# ╔═╡ b3caf772-2569-4ad3-a881-aaf5bf83bb4b
md"""
# Extracting sensitivity:

To extract the $\mu_U$ and calculate the sensitivty is as straight-forward as calculating the Credible interval of the posterior $p(n_S)$ and plugging the upper bound to the sensitivity equation: 

"""

# ╔═╡ d0959940-ff2a-439c-9352-5dd65f17762f
function get_tHalf_bayess( SNparams, eff, mu_U )
	return log(2)*SNparams["Nₐ"]*SNparams["foilMass"]*SNparams["a"]*SNparams["tYear"]/SNparams["W"]*eff/mu_U 
end;

# ╔═╡ b8c01575-a833-4010-9527-e6f2890be3c3
md"""
# Comparing the results: 

### Frequentist vs Bayess vs Cupid-0 (Bayess)
"""

# ╔═╡ 56e703f1-1b5b-4a96-a9f1-fcbc64e3746c
@bind clicked2 PlutoUI.Button("Generate data")

# ╔═╡ 1554ee0b-4a9c-49e4-badc-b000c767ed2d
data = let
	clicked
	clicked1
	clicked2
	generate_sample_data(h1d_full_bkg)
end

# ╔═╡ 23f06d08-e9ad-4de0-a152-fc5d47d6d103
h1d_sample = let 
	Hist1D(data; binedges=binningDict[:SumE]) 
end;

# ╔═╡ 26f30acb-b8ab-4619-bdb2-b49107c0a317
let
	f = Figure(size = (1050, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "sample background data \nlinear scale", xticklabelrotation=45)
	ax2 = Axis(f[1,2], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "log scale", yscale =log10, xticklabelrotation=45)
	
	st = hist!(ax, h1d_sample, label ="background")
	st2 = hist!(ax2, h1d_sample, label ="log background",)
	ylims!(ax2, 1e-1, 1e5)
	f
end

# ╔═╡ cff3ad43-3482-4705-9a85-f1b02558d641
md"""
!!! danger "Danger:" 
	The results using the general model vary **a lot**! A more robust model is necessary.
"""

# ╔═╡ b8715cfa-c2db-4905-9c06-ed37b784b98e
md"""
# Summary
### 1. A general comparison of Bayessian and frequentist statistics

- #### Key points: 
  - point-estimate ➡ pdf; 
  - confidence interval ➡ credible interval;
  - parameters fixed ➡ distributed

"""

# ╔═╡ bfd48054-c220-4cc3-b5ad-3f63a208a1eb
md"""
# Summary

### 2. Introduction to sensitivity
- #### Key points: 
  - "*how low the decay half-life must be in order to be able to measure it with our setup*"
  - Frequentist uses Feldman Cousins, depends only on $\bar{b}$
  - Bayessian calculates from posterior


"""

# ╔═╡ 47b130b9-9211-4756-a44c-1daedfb4e44a
md"""
# Summary
### 3. Simulated data-set
- #### Key points: 
  - Simulations using `falaise v5` and `simRC` 
  - Simulated a lot more events ➡ more precise even in high energy regions
  - Expected a lot more background than previously estimated ➡ problem?


"""

# ╔═╡ 4603c776-110e-4b48-82d5-314b61038fd7
md"""
# Summary
### 4. Model
- #### Key points: 
  - 2 possibilities for building likelihoods
  - **general model** gives very large $\bar{b}$ ➡ problem!
  - **detailed model** ➡ need to try out


"""

# ╔═╡ 5d1ab821-6b0b-42c2-a648-a0df740c114d
md"""
# Summary
### 5. Frequentist approach: "The Feldman Cousins approach"
- #### Key points: 
  - undesired sensitivities obtained ➡ maybe using better cuts (and tracking algo) will help?


"""

# ╔═╡ 783d55e2-82bb-4209-a349-2c0a001ab463
md"""
# Summary
### 6. Bayessian approach: "Here's what the data tells me..."
- #### Key points: 
  - depends on choice of priors (tried to make informative...but maybe I shouldn't?)
  - gives more optimistic results 🔥
  - actually can get close to cupid-0



"""

# ╔═╡ 0dc66ff5-c2c5-4b86-a971-f9f01e7dc7ff
md"""
# The end! 😏
"""

# ╔═╡ 66d151e6-09a3-4704-9e31-aaac862ab24b
md"""
# Backup
"""

# ╔═╡ 5bde50cd-62fe-4322-9826-61910bae0d04
md"""
# Using BAT.jl
"""

# ╔═╡ a83271ba-5664-4cd1-8058-7ef7f93a1e30
function likelihood_a(nB::Real, nS::Real, lambda::Real, x::Real; Q=2997.0, sigma = 60.0)
    (nB * pdf(Exponential(lambda), x)  +
     nS * pdf(Normal(Q, sigma),  x)) *
    inv( nB + nS ) 
end

# ╔═╡ 8b83a8d9-b723-40a1-8103-ee569534fb98

function likelihood_a(p::NamedTuple{(:nB, :nS, :lambda)}, x::Real; Q=2997.0, sigma = 60.0)
    likelihood_a(p.nB, p.nS, p.lambda, x; Q=Q, sigma = sigma)
end

# ╔═╡ aa951b6d-2d84-45f2-998c-92a8b062ac81
function likelihood_b(p, n)
    pdf(Poisson(p.nB + p.nS), n)
end

# ╔═╡ d1087de1-2533-425e-804d-800ca2f4ec87
struct UnbinnedModel
    data
end

# ╔═╡ b2eeec92-f2b5-469b-96f9-4dceba08650a
@inline DensityInterface.DensityKind(::UnbinnedModel) = IsDensity()

# ╔═╡ 4a2163e9-7931-4aea-9dd0-dba2d363a380
function DensityInterface.logdensityof(likelihood::UnbinnedModel, p)
    n = length(likelihood.data)
    nB, nS, lambda = p

    ll_a = sum(log.(likelihood_a.(nB, nS, lambda, likelihood.data)))
    ll_b = log(likelihood_b(p, n))
    ll_a + ll_b
end

# ╔═╡ 0c1f4549-0f10-4928-9864-daee87006f4c
prior = distprod(
    nB=Uniform(1e-5, 1e4), #Uniform(1e-5, 45),
    nS=Uniform(1e-15, 10), #Uniform(1e-5, 15),
    lambda = Uniform(1e-5,1e3)
)

# ╔═╡ ec996978-ea8d-47e6-b41c-3942e8105b23
function sci(
    X::AbstractVector{<:Real},
    W::AbstractWeights = UnitWeights{eltype(X)}(length(eachindex(X)));
    nsigma_equivalent::Real = 1
)
    nsigma_90percent = quantile(Normal(), 0.5 + 0.9/2)   # 90% = 1.6448536269514717

    m, n = if nsigma_equivalent ≈ oftype(nsigma_equivalent, 1)
        28,  41  # 0.6827 ≈ 28//41
    elseif nsigma_equivalent ≈ oftype(nsigma_equivalent, 2)
        42,  44  # 0.9545 ≈ 42//44
    elseif nsigma_equivalent ≈ oftype(nsigma_equivalent, 3)
        369,  370  # 0.9973 ≈ 369/370
    elseif isapprox(nsigma_equivalent, nsigma_90percent, atol = 0.01)   # 0.90 ≈ 1.64
        90, 100
    else
        throw(ArgumentError("nsigma_equivalent must be 1, 2, 3 or 1.64 (for 90% credibility interval)"))
    end

    qs = quantile(X, W, range(0, 1, length = n + 1))
    ivs = ClosedInterval.(qs[begin:end-1], qs[begin+1:end])
    
    sel_idxs = sort(sortperm(ivs, by = IntervalSets.width)[begin:begin+m-1])

    r_idxs = eachindex(sel_idxs)
    for i in r_idxs
        s = sel_idxs[i]
        if (i == first(r_idxs) || sel_idxs[i-1] != s-1) && (i == last(r_idxs) || sel_idxs[i+1] != s+1)
            if i >= first(r_idxs)+2 && sel_idxs[i-1] == s-2 && IntervalSets.width(ivs[s-1]) >= IntervalSets.width(ivs[s])/2
                sel_idxs[i] = s-1
            end
            if i <= last(r_idxs)-2 && sel_idxs[i+1] == s+2 && IntervalSets.width(ivs[s+1]) >= IntervalSets.width(ivs[s])/2
                sel_idxs[i] = s+1
            end
        end
    end

    sel_ranges = UnitRange{Int}[]
    for i in sel_idxs
        if isempty(sel_ranges) || i != last(sel_ranges[end]) + 1
            push!(sel_ranges, i:i)
        else
            sel_ranges[end] = first(sel_ranges[end]):i
        end
    end

    [ClosedInterval(minimum(ivs[first(r)]), maximum(ivs[last(r)])) for r in sel_ranges]
end

# ╔═╡ 21253d0e-0676-4615-82d4-13d52aad20c8
md"""
### Changing the prior to "more informative" (Normal around the MLE values)
"""

# ╔═╡ f7107b98-4e2e-418c-8f5e-773decc18eda
md"""
### Comparison with FWHM = 12%
"""

# ╔═╡ e2161672-3e61-4163-a7b6-5cf0a075462c
allProcesses_12 = load_processes("fal5_12pers")

# ╔═╡ 6c93c4a9-8c02-4de2-adf1-02d2c642c8bb
begin
	sig12_process = deepcopy(get_process("bb0nu_foil_bulk", allProcesses_12))
	set_nTotalSim!(sig12_process, 1e8)
	nothing
end

# ╔═╡ 3483d64d-0e77-44d1-8819-e508daf5ae78
begin
	sig_process_12 = deepcopy(get_process("bb0nu_foil_bulk", allProcesses_12))
	set_nTotalSim!(sig_process_12, 1e8)
	nothing
end

# ╔═╡ 44912c3b-7bf9-4245-9f8e-32725dec45a5
begin
	bb2nuProcess_12 = deepcopy(get_process("bb_foil_bulk", allProcesses_12))
	set_signal!(bb2nuProcess_12, false)
	bkgProcessesESum_12 = [
	    bb2nuProcess_12,
	    deepcopy(get_process("Bi214_foil_bulk", allProcesses_12)),
	    deepcopy(get_process("Bi214_foil_surface", allProcesses_12)),
	    deepcopy(get_process("Bi214_wire_surface", allProcesses_12)),
	    deepcopy(get_process("Tl208_foil_bulk", allProcesses_12)),
	]
	set_nTotalSim!.(bkgProcessesESum_12, 1e8)
	
	nothing
end

# ╔═╡ 348c76b8-6207-415d-97e5-69004a34da2f
expBkgESum_12, effbb_12 = let
	t12MapESum = get_tHalf_map(SNparams, α, sig_process_12, bkgProcessesESum_12...)
	best_t12ESum = get_max_bin(t12MapESum)
	expBkgESum = get_bkg_counts_ROI(best_t12ESum, bkgProcessesESum_12...)
	effbb = lookup(sig_process_12, best_t12ESum)
	ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

	expBkgESum, effbb
end

# ╔═╡ 3dc657f6-692b-47ed-8c4c-2d3eec90d95e
get_tHalf2(t, b, eff) = get_tHalf(
    SNparams["W"],
    SNparams["foilMass"],
    SNparams["Nₐ"],
    t,
    SNparams["a"],
    eff,
    b / SNparams["tYear"] * t,
    α;
    approximate="formula"
)

# ╔═╡ 78d57fb7-f166-420b-af5d-7f982252a0b4
md"""
# Auxiliary code blocks:
"""

# ╔═╡ 5d6fb41a-8146-4423-94b2-dc2e78495a25
f_exp(x, p) = p[1] * exp(-x * p[2]);

# ╔═╡ 3539584e-19b7-45d1-b152-2a0d45e03755
bw = step(binningDict[:SumE]);

# ╔═╡ e881e2dc-a9ac-4ea2-ab29-9da4d4cc8b2f
function fit_function(x, p; Q=2997.0-a, sigma = 60, deltaE=bw)
	lambda, nS, nB = p
    deltaE*(nB * pdf( Exponential(inv(lambda)), x ) + nS * pdf( Normal(Q, sigma), x) )
end

# ╔═╡ 095600d3-bf7d-4373-9e92-b92557010b49
# data needs to be shifted to start from 0 for a binned fit (we do this by subtracting a from edges)
sample_ROI_data = filter( x-> a <= x <= b, data ) .- a;

# ╔═╡ 9e0d7d49-a7d2-458a-b8f4-c3b69e7b9e8e
unbinned_posterior = PosteriorMeasure(UnbinnedModel(sample_ROI_data), prior)


# ╔═╡ 0e0841e6-62b9-484e-8eb9-07c2d4de4652
unbinned_samples = bat_sample(unbinned_posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result


# ╔═╡ 05c00d88-b2b4-4c6f-8e45-1842681b12f2
unbinned_unshaped_samples, f_flatten = bat_transform(Vector, unbinned_samples)


# ╔═╡ 85d7f83f-4fce-49bb-9df8-82a37aebac18
muSsUnbinned = [m[2] for m in unbinned_unshaped_samples.v]


# ╔═╡ aa299272-0ff1-497c-95bf-66c122594489
n_S_BAT = IntervalSets.rightendpoint(sci(muSsUnbinned, nsigma_equivalent=1.65)[1])

# ╔═╡ eddc431c-0560-4cc1-a5e5-93c0459eb9ea
lambda = fit_mle(Exponential, sample_ROI_data).θ |> inv;

# ╔═╡ 8ba6e3d9-3529-4c8d-bc5c-98a91743f29a
@model function full_model_exponential(data; Q=2997.9-a, sigma=60)
	lambda ~ Uniform(1e-3, 1e3)
    # Prior distributions for the proportions
    nSignal ~ Uniform(0.0, 10.0)  # Prior for proportion of Normal data
    nBkg ~ Uniform(0.0, 1e3)  # Prior for proportion of Normal data
    
    # Likelihood of data
    for x in data
        likelihood_normal = pdf(Normal(Q, sigma), x)
        likelihood_exponential = pdf(Exponential(inv(lambda)), x)
        likelihood_mixture = (nSignal * likelihood_normal + nBkg * likelihood_exponential) 
        likelihood_extended = pdf(Poisson(nSignal + nBkg), length(data))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end;

# ╔═╡ c145389c-b760-42d2-afe2-18ea3d31b49f
@model function full_model_exponential_1(data; Q=2997.9-a, sigma=60)
	lambda ~ truncated(Normal(1/100.0, 1.0), 0.0, Inf)
    # Prior distributions for the proportions
    nSignal ~ Uniform(0.0, 10.0)  # Prior for proportion of Normal data
    nBkg ~ Uniform(0.0, 1e3)  # Prior for proportion of Normal data
    
    # Likelihood of data
    for x in data
        likelihood_normal = pdf(Normal(Q, sigma), x)
        likelihood_exponential = pdf(Exponential(inv(lambda)), x)
        likelihood_mixture = (nSignal * likelihood_normal + nBkg * likelihood_exponential) 
        likelihood_extended = pdf(Poisson(nSignal + nBkg), length(data))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end;

# ╔═╡ e4822025-24e7-4390-a98f-6b2137174d2c
chains2 = Turing.sample(full_model_exponential_1(sample_ROI_data), NUTS(0.65), MCMCThreads(), 5_000, 2 );

# ╔═╡ 33fe2df3-81e5-4b1c-8930-56651e09228d
N = bw * length(sample_ROI_data) * lambda;

# ╔═╡ 884b29e7-9cf3-4779-9b04-ffc35ebaff55
fit_params = [N, lambda];

# ╔═╡ c48ff5b5-cfcc-48b1-9859-70d5c0fae260


# ╔═╡ a8025e0c-a17a-467f-9873-6999847d184f
h1d_sample_ROI = restrict(h1d_sample, a, b);

# ╔═╡ e20a262d-362b-40d1-83a3-826525251f50
function chi_square(bincounts, f, xs; nparams= 1)
	ndf = length(bincounts) - nparams
	χ² = 0.0
	for i in eachindex(bincounts)
		if ( bincounts[i] >= 5 )
			χ² += (f(xs[i] - xs[1]) - bincounts[i])^2 / f(xs[i] - xs[1])
		else
			ndf -= 1
		end
	end
	
	return χ², ndf
end;

# ╔═╡ 8f81bf60-cb6a-4d24-a793-18cb621c1afb
let l = lambda, N = N, h = h1d_sample_ROI, fit_f(d) = f_exp(d, fit_params)
	
	f = Figure(size = (800, 500), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel ="energy", ylabel="counts", xticks= range(a, b, 5), yscale =log10)
	
	bin_edges = binedges(h)
	n = sum(bincounts(h))
	min_bin, max_bin = extrema(bin_edges) # used for scaling on the x-axis
	
	xs = range(min_bin, max_bin, length=100) 
	ys = fit_f.(xs .- min_bin) 

	chisq, ndf = chi_square( bincounts(h), fit_f, bincenters(h))

	p = plot!(ax, h, label = "data", xlabel ="energy", ylabel="normalized counts", )
	vspan!(ax, best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge], color = (:orange, 0.3), label ="ROI" )
	lines!(ax, xs , ys, color = :red, linewidth=3, linestyle=:solid, label ="fit: \nλ = $(round(l, sigdigits=3)) \nχ² / ndf = $(round(chisq, sigdigits= 4)) / $(ndf)")
	errorbars!(ax, h, label="data", whiskerwidth = 6, color=:black, )
	ylims!(ax, 1e-1, 1e4)
	f[1,2] = Legend(f, ax, merge = true )
	f
end

# ╔═╡ 79b85562-16a8-423d-ac7d-16d9ab89c61b
function exp_integral(lambda, N, a, b)
	first = -inv(lambda)*N*exp(-a*lambda)
	second = -inv(lambda)*N*exp(-b*lambda)
	
	return second - first
end;

# ╔═╡ 890ef7fc-2649-442e-beaf-8f9db8c69220
# data needs to be shifted to start from 0 for a binned fit (we do this by subtracting a from edges)
b_exp_sample = exp_integral(lambda, N, best_t12ESum[:minBinEdge] - a  , best_t12ESum[:maxBinEdge] - a) / bw ;

# ╔═╡ df989bf7-ee2c-4912-bee8-42afeeca70e1
T12_sample_freq = round(get_tHalf(SNparams, effbb, b_exp_sample, α), sigdigits = 4);

# ╔═╡ ef493c16-656c-4253-bb46-62e0aefae4c9
md"""
# From this fit we can extract expected number of background counts in ROI per 17.5kgy. 

### Estimated $\bar{b} =$ $(round(b_exp_sample, sigdigits=4)) using the **general model**.

#### This value can be plugged into the sensitivity equation to obtain: $T^{1/2} \geq$ $T12_sample_freq y

!!! danger "Note on the result:" 
	- The result of the fit changes drastically with the choice of fit region and fluctuation of data itself (generating different sample data-set changes results!)
	- This is why a **proper** background model and choice of likelihood is necessary! 
	- Fitting with an exponential is not ideal, the **detailed model** should be used instead. 
"""

# ╔═╡ dc1f97b4-8e33-42ca-a434-3a6cdc22abf9


# ╔═╡ ee169e16-4cc7-4a3b-b5cd-58a4fc0bfe2f
function likelihood( params, x )
	Q, sigma, lambda, nSig, nBkg = params
	if( lambda <= 0.0 )
		return 0.0
	else
    	return (nSig * pdf( Normal(Q, sigma), x ) + nBkg * pdf(Exponential(lambda), x)) / ( nSig+nBkg )
	end
end;

# ╔═╡ eea1764c-e8f1-4fe2-ba00-a497ffbe3462
function Makie.plot(ch::Chains)
	fig = Figure(size = (1200, 400), fontsize= 22, fonts = (; regular = "TeX"))
	for (ind, param) in enumerate(ch.name_map.parameters)
		ax = Axis(fig[1, ind], title=string(param), xticklabelrotation = 45, xlabel ="sample value", ylabel = "counts" )
		for (ind2, datavec) in enumerate(eachcol(getindex(ch, param).data))
			# Get current default colorpalette
			colors = Makie.current_default_theme().attributes[:palette][][:color][]
			stephist!(ax, datavec, color=(colors[ind2], 1.0), bins = 30,
				linewidth = 2,
			)
			xlims!(ax, extrema(datavec))
		end
	end
	return fig
end

# ╔═╡ 5261a274-7db2-47c1-bf0b-ba0bcaf2727a
plot(chains2)

# ╔═╡ 6ad5eb74-73d9-4ee4-8383-e47a783bc6f9
model = full_model_exponential(sample_ROI_data);

# ╔═╡ 05af0564-6522-41ff-80e5-41753ce572ed
begin
	train
	chains = Turing.sample(model, NUTS(0.65), MCMCThreads(), 10_000, 2, burn_in=1000, thinning=2 );
	nothing
end

# ╔═╡ 907872b9-4aa7-4c52-9473-d4d2a4b77fc1
summarystats(chains)

# ╔═╡ 1a6c1f12-17fb-448d-86ec-bfa2906a34d7
let
	f = plot(chains)
	f
end

# ╔═╡ fee38023-7b43-42c3-b74b-da604d72dfba
l68_l, l68_u  = round.(quantile(chains[:lambda], [0.5-0.68/2, 0.5+0.68/2]), sigdigits = 3);

# ╔═╡ 0b121bf6-b573-4b05-9387-bd03b5d8f9e9
md"""
### Here we see the 3 posterior distributions, after sampling, of the parameters of interest. Notice:
-  $p(\lambda)$ has the peak close to MLE fit value of $\lambda_{MLE} =$ $(round(lambda, sigdigits= 3)); the uncertainty can be evaluated from the posterior distribution
  - i.e. 68% CI ∈ ( $l68_l, $l68_u )
-  $p(n_B)$ has the peak close to the total number of events in the sample set: $(length(sample_ROI_data))
-  $p(n_S)$ follows a roughly **exponential** shape with mode at 0, i.e. 0 signal counts. **From this distribution we estimate $\mu_U$!**
"""

# ╔═╡ c65c05b9-835f-454d-a618-b3c9b75cf537
n_S_turing = IntervalSets.rightendpoint(sci(Array(chains[:,2,1]), nsigma_equivalent=1.65)[1])

# ╔═╡ f07e5dbf-84f0-4ae9-8437-6e1c9bda29ee
t12_bayess_2 = let 
	n_Sig_90q = MCMCChains._hpd(Array(chains[:, 2, 1]); alpha=0.1)[2]
	t12_bayess = get_tHalf_bayess(SNparams, effbb, n_Sig_90q)
	t12_bayess
end

# ╔═╡ e00492d8-872d-4e0c-a067-96e0c885c28f
begin
	λ = mode(chains[:,1,:])
	n_Sig = mode(chains[:,2, :]) 
	n_Sig_90q = MCMCChains._hpd(Array(chains[:, 2, 1]); alpha=0.1)[2]
	
	n_Bkg = mode(chains[:,3, :]) 
	n_Bkg_90q = quantile(chains[:, 3, 1], 0.9);
	nothing
end

# ╔═╡ 948d9a83-a0a8-450e-95aa-a4c7775402a5
let l = λ, nS = n_Sig, nB = n_Bkg, h = h1d_sample_ROI, ch = chains

	f = Figure(size = (900, 700), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel ="energy [keV]", ylabel="counts / 100keV", title= "Bayessian fit \nlog-scale",  xticks= range(a, b, 5), yscale = log10)
	ax2 = Axis(f[2,1], xlabel ="energy [keV]", ylabel="counts / 100keV", title= "Bayessian fit \nlinear-scale",  xticks= range(a, b, 5), )
	p = plot!(ax, h, label = "data")
	p = plot!(ax2, h, label = "data")

	bin_edges = binedges(h)
	min_bin, max_bin = extrema(bin_edges) # used for scaling on the x-axis
	
	xs = range(min_bin, max_bin, length=1000) 
	colors = [:red, :green, :yellow]
	alphas=[0.01, 0.10, 1-0.68]
	m, l, u = get_mean_conf(xs .- min_bin, ch, fit_function, alphas[1])	
	for i=1:3
		m, l, u = get_mean_conf(xs .- min_bin, ch, fit_function, alphas[i])	
		band!(ax, xs, m .-l, u .+m ; color = (colors[i], 1), label="$(1-alphas[i])% CI")
		band!(ax2, xs, m .-l, u .+m ; color = (colors[i], 1), label="$(1-alphas[i])% CI")
	end
	med = fit_function( xs .- min_bin, [λ,n_Sig,n_Bkg]; Q = 2997.3 )
	lines!(ax, xs, med, color = (:black), label="Bayess median")
	lines!(ax2, xs, med, color = (:black), label="Bayess median")
	errorbars!(ax, h, label="data", whiskerwidth = 6, color=:black, )
	errorbars!(ax2, h, label="data", whiskerwidth = 6, color=:black, )
	ylims!(ax, 1e-2, 1e3)
	f[1:2,2] = Legend(f, ax, merge = true )
	colgap!(f.layout, 1, Relative(0.07))
	f


end

# ╔═╡ 0ba2b76d-40c9-4846-8f8b-aeb5ca4a9b33
t = range(0, 6, 60);

# ╔═╡ e5328719-4e51-4968-a124-fc49d61bac4d
get_tHalf1(t, b) = get_tHalf(
    SNparams["W"],
    SNparams["foilMass"],
    SNparams["Nₐ"],
    t,
    SNparams["a"],
    effbb,
    b / SNparams["tYear"] * t,
    α;
    approximate="formula"
);

# ╔═╡ 6b91b99c-c299-4174-afbd-9dcd85ff0a27
t12_bayess = log(2)*SNparams["Nₐ"]*SNparams["foilMass"]*SNparams["a"]*SNparams["tYear"]/SNparams["W"]*effbb/n_Sig_90q ;

# ╔═╡ af325d7d-6c1b-475b-b608-8300d8fc78e5
let 
	fig = Figure( size = (1200, 500), fontsize= 22, fonts = (; regular = "TeX") )
	ax1 = Axis(fig[1,1] , limits = (0, 10, nothing, nothing), xlabel =L"\textrm{sample value}", ylabel = L"\textrm{pdf}", title = L"p(n_S|data)\textrm{ posterior distribution}")
	ax2 = Axis(fig[1,2], yscale = log10, xlabel =L"\textrm{\mu_U}", ylabel = L"\textrm{sensitivity [y]}", title =L"\textrm{sensitivity as a function of upper bound on n_S}")

	h = normalize(Hist1D(chains[:,2, 1]; binedges = 0:0.1:10))
	h1 = restrict(h, 0, n_Sig_90q)
	xs = bincenters(h1)
	ys = bincounts(h1)
	c1 = Makie.wong_colors()[5]
	p1 = stairs!(ax1, h, color = c1, label = L"\textrm{posterior samples}")
	hist!( ax1, h1, color = c1, strokewidth = 0, strokecolor = c1, label = L"\textrm{90% credible interval}" )
	axislegend(ax1, merge = true, position=:rt, backgroundcolor= :white, margin= (0,15,10,15))
	ax1.xticks= ([0, n_Sig_90q, 5, 10], ["0", string(round(n_Sig_90q, digits =2)), "5", "10"])	

	n_Ss = sort(chains[:,2, 1])
	t12s = map( x -> get_tHalf_bayess( SNparams, effbb, x), n_Ss) 
	t12s_excl = filter( x -> x <=  t12_bayess, t12s )  
	n_excl = length(t12s_excl)
	lines!(ax2, n_Ss,t12s)
	n_Ss[end-length(t12s_excl)+1:end][1]
	band!( ax2, n_Ss[end-n_excl+1:end], 1e23 .*ones(length(t12s_excl)), t12s_excl, color =(:red, 0.4) )
	band!( ax2, n_Ss[1:Int(length(t12s)-n_excl)+1], 1e23 .*ones(Int(length(t12s)-n_excl)+1),  t12_bayess.*ones(Int(length(t12s)-n_excl)+1), color =(:red, 0.4) )
	text!( ax2, 2.5, 5e23, text=("excluded region") )
	ylims!(ax2,1e23, 1e28)
	xlims!(ax2,0, 8)
	colsize!(fig.layout, 1, Relative(0.4))
	fig
end

# ╔═╡ 6599d650-f28e-416d-aeb1-6afb0e7c3672
md"""
### The results are:	
####  ``\mu_U`` = $(round(n_Sig_90q, sigdigits =3)) 
####  ``T^{1/2} \geq`` $(round(t12_bayess ,sigdigits=3)) yr

!!! update " " 
	The obtained estimate is better than the one obtained from frequentist approach! 
"""

# ╔═╡ f45babc7-ba7f-42f3-9192-521ee77058d7
let 
	fig = Figure( size = (1100, 500), fontsize= 22, fonts = (; regular = "TeX") , italic="TeX")
	ax1 = Axis(fig[1,1] , limits = (0, 6, nothing, nothing), ylabel =L"\textrm{sensitivity [yr]}", xlabel = L"\textrm{measurement duration [y]}", title = L"\textrm{sensitivity at 90% CL for SN and Cupid-0}")
	
	hlines!(ax1, 4.6e24, color = :black, label = L"\textrm{Cupid: 4.6\times10^{24}y}", linestyle =:dash, linewidth=4)
	lines!(ax1, t, get_tHalf1.(t, expBkgESum), label =L"\textrm{freq: 8% simulation; } \bar{b} = %$(round(expBkgESum, sigdigits=3))", linewidth=4)
	lines!(ax1, t, get_tHalf1.(t, b_exp_sample), label =L"\textrm{freq: 8% sample set; } \bar{b} = %$(round(b_exp_sample, sigdigits=3))", linewidth=4)
	#lines!(ax1, t, get_tHalf2.(t, expBkgESum_12, effbb_12), label =L"\textrm{freq: 12% baseline}", linewidth=4, linestyle =:dot)
	scatter!(ax1, [2.88], [t12_bayess], color = :red,label = L"\textrm{Bayess: 8% sample set; } \mu_U = %$(round(n_Sig_90q, sigdigits=3))", marker= :star5, markersize = 25)
	ylims!(ax1, 0, 6e24)
	fig[1,2] = Legend(fig, ax1, merge = true, tellwidth = false,patchsize = (40, 30))
	colgap!(fig.layout, 1, Fixed(-10))
	fig
end

	

# ╔═╡ dfa0066b-1040-4575-9cde-d040d8e3cdce
md"""
#### with informative prior: T12 = $t12_bayess
"""

# ╔═╡ ba3142e4-f68e-48a1-9dfd-2fa4695a1488
let 
	fig = Figure( size = (1000, 500), fontsize= 22, fonts = (; regular = "TeX") )
	ax1 = Axis(fig[1,1] , limits = (0, 6, nothing, nothing), xlabel =L"\textrm{sensitivity [yr]}", ylabel = L"\textrm{measurement time}", title = L"\textrm{sensitivity at 90% CL for SN and Cupid-0}")
	
	lines!(ax1, t, get_tHalf1.(t, expBkgESum), label =L"\textrm{freq: 8% baseline}", linewidth=4)
	lines!(ax1, t, get_tHalf2.(t, expBkgESum_12, effbb_12), label =L"\textrm{freq: 12% baseline}", linewidth=4, linestyle =:dot)
	lines!(ax1, t, get_tHalf1.(t, b_exp_sample), label =L"\textrm{freq: 8% general model}", linewidth=4)
	scatter!(ax1, [2.88], [t12_bayess], color = :red,label = L"\textrm{Bayess: 8% general model}", marker= :star5, markersize = 25)
	hlines!(ax1, 4.6e24, color = :black, label = L"\textrm{Cupid: 4.6\times10^{24}y}", linestyle =:dash, linewidth=4)
	fig[1,2] = Legend(fig, ax1, merge = true )
	fig
end

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
# ╟─503673d7-7dac-43eb-a493-0aea57674720
# ╟─50e7af06-2bce-4914-a867-ec6261a6714b
# ╟─839d78ac-2ad9-4cc0-92df-7d3402a5b291
# ╟─29274571-9856-4dd7-9afa-d5204fac46ca
# ╟─50b657c7-f326-4cfa-a4dc-160980d1ec24
# ╟─80b557d9-96b8-4d00-93cd-2986798a6054
# ╟─93d6b82b-6d0e-4be5-812b-5cee189987f2
# ╟─7a925206-47b3-4702-9df2-4074239ca4e6
# ╟─ce031aaf-ceb5-4b82-8e0e-bc3a795daecd
# ╟─50493bc0-eb22-4af5-ac47-5da71eb96632
# ╟─74d7f5c7-d5bb-41c1-9e72-b5da45cb3538
# ╟─cbf182c3-4e2e-445d-9814-f5f6146a25e4
# ╟─7f8276ae-0a48-40ad-bdec-ee2cefc7ec55
# ╟─7be26605-60c2-4a03-a090-385caa299091
# ╟─68da87d9-c5bb-4aa5-837b-6803fb7687e6
# ╟─ee33bb4c-3f6b-4be3-be77-18e121208444
# ╟─3a34adbb-bc3c-4043-a524-e0697523900d
# ╟─bee405ae-e2bd-48b8-8de8-652b99bbd8c4
# ╟─955810a1-9ef6-4531-8156-02b070c29676
# ╟─11b39aa3-2c3b-4d5f-8ef0-e938c3cb3399
# ╟─e4c6c2b0-2afe-4667-894c-71c8476239c6
# ╟─94ad3e49-34d8-45e5-ac65-a9241136f226
# ╟─875130d7-f0e3-4be0-a963-c427712ccdcf
# ╟─557e5713-6563-4e4e-b0ce-a9ee67cf12e6
# ╟─b1991653-1363-4810-a78d-923c4ce64daa
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
# ╟─68231ef5-9425-4e4c-98a8-586323938d45
# ╟─1c44a17c-41a5-4ec2-a37e-b4c20e9407cc
# ╟─0033c567-5cf4-4067-bc50-00543dabc798
# ╟─c69bf1d5-e7ae-4b31-8fe3-cbb872e68424
# ╟─58241414-5ec2-463a-8422-f7bd8e30b606
# ╟─91546e37-2d22-4616-9d19-a4d973e80e54
# ╟─d5b67056-4212-40b6-b25c-007cacd0509a
# ╟─26f30acb-b8ab-4619-bdb2-b49107c0a317
# ╟─06e1024a-5afb-4cd8-a75e-4bcfb4a5ab8f
# ╟─1554ee0b-4a9c-49e4-badc-b000c767ed2d
# ╟─23f06d08-e9ad-4de0-a152-fc5d47d6d103
# ╟─a47c9336-33de-4c89-bc64-3e1e12036243
# ╟─873c693a-23a9-4ec2-ac3a-04cd9c9ac40b
# ╟─52e4668e-9804-492a-9deb-830f8cdedfcf
# ╟─fab2761b-fa89-40e1-bfd8-e30d5524b7c2
# ╟─f768a9e3-7ab7-4c36-bf5b-2737a638c325
# ╟─370440df-155a-46c2-a3f1-75f4ebf73bc9
# ╟─255a4f62-964b-47d6-96f0-663cfcdcf09f
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
# ╟─6ea1a763-3268-4f77-83a0-f6065d77eb0a
# ╟─69f34564-f7e5-45d5-90d6-501591ae67ae
# ╟─2485d73a-28bc-461a-be00-0083dba8a5fb
# ╟─796d7a86-dae3-4e0c-9dc7-157d245425e9
# ╟─8f81bf60-cb6a-4d24-a793-18cb621c1afb
# ╟─ef493c16-656c-4253-bb46-62e0aefae4c9
# ╟─46086903-ec96-4f86-8081-ddfb713fde30
# ╟─d941baf4-d60a-40f8-a35d-bf76e16dbeda
# ╟─742f7fd5-01e7-4731-a7a1-ad84bb86b8c7
# ╟─278988c0-6da6-46f2-a50d-89ea32edb69a
# ╟─4cd71a16-3bb8-4d81-b1f7-f67eb6cad196
# ╟─14ef85c1-4be6-45c3-b2f8-cae5d5a21aa7
# ╟─2c04e398-8c31-4f9b-a00a-528c9bf60358
# ╟─8ba6e3d9-3529-4c8d-bc5c-98a91743f29a
# ╟─c145389c-b760-42d2-afe2-18ea3d31b49f
# ╟─6376b61a-f6a2-4dbd-9ba1-b2bcd96f5d9d
# ╟─2134a77e-34bb-4b7f-a268-23fd1e6d52a1
# ╟─b9047089-49d5-4ad7-8e2c-d2029176b393
# ╟─05af0564-6522-41ff-80e5-41753ce572ed
# ╟─7ca9d29b-688b-45a3-bb13-dce16e067fae
# ╟─907872b9-4aa7-4c52-9473-d4d2a4b77fc1
# ╟─14bde55d-7695-431c-930f-2b7181665937
# ╟─0b121bf6-b573-4b05-9387-bd03b5d8f9e9
# ╟─1a6c1f12-17fb-448d-86ec-bfa2906a34d7
# ╟─fee38023-7b43-42c3-b74b-da604d72dfba
# ╟─dcc61c0b-5529-42b7-8fd1-88033b060b0b
# ╟─9a8cff73-af96-4c5b-a455-4b890f9c3dc3
# ╟─948d9a83-a0a8-450e-95aa-a4c7775402a5
# ╟─c2bdf038-475e-4b52-8c4a-92e1ad09c1e0
# ╟─b3caf772-2569-4ad3-a881-aaf5bf83bb4b
# ╟─af325d7d-6c1b-475b-b608-8300d8fc78e5
# ╟─6599d650-f28e-416d-aeb1-6afb0e7c3672
# ╟─d0959940-ff2a-439c-9352-5dd65f17762f
# ╟─b8c01575-a833-4010-9527-e6f2890be3c3
# ╟─f45babc7-ba7f-42f3-9192-521ee77058d7
# ╟─56e703f1-1b5b-4a96-a9f1-fcbc64e3746c
# ╟─cff3ad43-3482-4705-9a85-f1b02558d641
# ╟─b8715cfa-c2db-4905-9c06-ed37b784b98e
# ╟─bfd48054-c220-4cc3-b5ad-3f63a208a1eb
# ╟─47b130b9-9211-4756-a44c-1daedfb4e44a
# ╟─4603c776-110e-4b48-82d5-314b61038fd7
# ╟─5d1ab821-6b0b-42c2-a648-a0df740c114d
# ╟─783d55e2-82bb-4209-a349-2c0a001ab463
# ╟─0dc66ff5-c2c5-4b86-a971-f9f01e7dc7ff
# ╟─66d151e6-09a3-4704-9e31-aaac862ab24b
# ╠═5bde50cd-62fe-4322-9826-61910bae0d04
# ╠═d33cdccb-3c3e-4394-b16b-6a0bd56b5a3f
# ╠═a83271ba-5664-4cd1-8058-7ef7f93a1e30
# ╠═8b83a8d9-b723-40a1-8103-ee569534fb98
# ╠═aa951b6d-2d84-45f2-998c-92a8b062ac81
# ╠═d1087de1-2533-425e-804d-800ca2f4ec87
# ╠═b2eeec92-f2b5-469b-96f9-4dceba08650a
# ╠═4a2163e9-7931-4aea-9dd0-dba2d363a380
# ╠═0c1f4549-0f10-4928-9864-daee87006f4c
# ╠═9e0d7d49-a7d2-458a-b8f4-c3b69e7b9e8e
# ╠═0e0841e6-62b9-484e-8eb9-07c2d4de4652
# ╠═05c00d88-b2b4-4c6f-8e45-1842681b12f2
# ╠═85d7f83f-4fce-49bb-9df8-82a37aebac18
# ╠═aa299272-0ff1-497c-95bf-66c122594489
# ╠═c65c05b9-835f-454d-a618-b3c9b75cf537
# ╠═ec996978-ea8d-47e6-b41c-3942e8105b23
# ╟─21253d0e-0676-4615-82d4-13d52aad20c8
# ╟─e4822025-24e7-4390-a98f-6b2137174d2c
# ╟─5261a274-7db2-47c1-bf0b-ba0bcaf2727a
# ╠═f07e5dbf-84f0-4ae9-8437-6e1c9bda29ee
# ╟─dfa0066b-1040-4575-9cde-d040d8e3cdce
# ╟─f7107b98-4e2e-418c-8f5e-773decc18eda
# ╠═e2161672-3e61-4163-a7b6-5cf0a075462c
# ╠═6c93c4a9-8c02-4de2-adf1-02d2c642c8bb
# ╠═3483d64d-0e77-44d1-8819-e508daf5ae78
# ╠═44912c3b-7bf9-4245-9f8e-32725dec45a5
# ╠═348c76b8-6207-415d-97e5-69004a34da2f
# ╠═ba3142e4-f68e-48a1-9dfd-2fa4695a1488
# ╠═3dc657f6-692b-47ed-8c4c-2d3eec90d95e
# ╟─78d57fb7-f166-420b-af5d-7f982252a0b4
# ╠═e881e2dc-a9ac-4ea2-ab29-9da4d4cc8b2f
# ╠═a5b84247-8464-4a80-af07-5309d7bcaef9
# ╠═5d6fb41a-8146-4423-94b2-dc2e78495a25
# ╠═3539584e-19b7-45d1-b152-2a0d45e03755
# ╠═eddc431c-0560-4cc1-a5e5-93c0459eb9ea
# ╠═33fe2df3-81e5-4b1c-8930-56651e09228d
# ╠═884b29e7-9cf3-4779-9b04-ffc35ebaff55
# ╠═df989bf7-ee2c-4912-bee8-42afeeca70e1
# ╠═890ef7fc-2649-442e-beaf-8f9db8c69220
# ╠═095600d3-bf7d-4373-9e92-b92557010b49
# ╠═c48ff5b5-cfcc-48b1-9859-70d5c0fae260
# ╠═a8025e0c-a17a-467f-9873-6999847d184f
# ╠═e20a262d-362b-40d1-83a3-826525251f50
# ╠═79b85562-16a8-423d-ac7d-16d9ab89c61b
# ╠═dc1f97b4-8e33-42ca-a434-3a6cdc22abf9
# ╠═ee169e16-4cc7-4a3b-b5cd-58a4fc0bfe2f
# ╠═eea1764c-e8f1-4fe2-ba00-a497ffbe3462
# ╠═6ad5eb74-73d9-4ee4-8383-e47a783bc6f9
# ╠═e00492d8-872d-4e0c-a067-96e0c885c28f
# ╠═752369c9-464f-44fc-97df-0f642bc83af7
# ╠═0ba2b76d-40c9-4846-8f8b-aeb5ca4a9b33
# ╠═e5328719-4e51-4968-a124-fc49d61bac4d
# ╠═6b91b99c-c299-4174-afbd-9dcd85ff0a27
