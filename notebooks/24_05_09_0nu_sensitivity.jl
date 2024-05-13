### A Pluto.jl notebook ###
# v0.19.40

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
	include(scriptsdir("LoadData.jl"))
end

# ╔═╡ 039334ac-d249-4c5f-bc97-e1180b52b47b
using LsqFit

# ╔═╡ a5b84247-8464-4a80-af07-5309d7bcaef9
using PlutoUI

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
### SuperNEMO Analysis meeting 9.5.2024

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


# ╔═╡ bee405ae-e2bd-48b8-8de8-652b99bbd8c4
md"""
# 2. Introduction to sensitivity
"""

# ╔═╡ 955810a1-9ef6-4531-8156-02b070c29676
md"""
# But what is sensitivity really? 

Imagine an ideal experiment with 0 background whatsoever. In this "simple" experiment we are trying to estimate the half-life of a rare decay. To do so we can use the following forumla:

```math 
T^{1/2} = ln(2)\frac{N_A \cdot m \cdot t }{W}\frac{\varepsilon}{n_{S}} 
```

!!! danger "Issues with this equation:"
	- If the process is allowed (remember, 0nu could really be forbidden) we expect to measure some number of events $n$. 
	-  $n$ does not necessarily have to equal the **expected number of signal events** $n_S$.
	- Radioactive decay is a statistical process - we may be lucky and get many events, or we may be unlucky and get nothin'...
	- What do we do if the number of observed **counts** is $n = 0$? Infinite half-life? 🤔
"""

# ╔═╡ 11b39aa3-2c3b-4d5f-8ef0-e938c3cb3399
md"""
### Statistics to the rescue! (It is both our friend and our enemy...)

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
# Now about the $\mu_U$... 

#### There are two ways of looking at what it means:
1. ## Bayessian:
!!! note "Bayessian terminology:" 
	- In Bayessian analysis, $n_S$ is a parameter of intereset and is distributed according to some pdf $p(n_S | data)$. 
	- A **credible interbal** is a range of plausible values that $n_S$ can take, i.e. a 90% credible interval covers 90% of plausible values. (The choice of central, one-sided, highest density is up to the scientist.) 
	- 
	 **$\mu_U$ is then the upper bound of the credible interval!**

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
	ax.xticks= ([0, 2.7], ["0", L"\mu_U"])
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

### Here, I will show an example of calculating sensitivity to $0\nu\beta\beta$ for SuperNEMO using both approaches, Bayessian and frequentist. I will do so on a generated *sample data-set*. 

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
7. Background: `Se82.2nubb`, `Bi214_foil_bulk`, `Bi214_foil_surface`, `Bi214_wire_surface`, `Tl208_foil_bulk`


"""

# ╔═╡ cede3473-fae8-4998-ba6e-01f9d05ad54c
md"""
# Here's what the raw data looks like:
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
	ax = Axis(f[1:2,1], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "stacked normalized histogram", xticklabelrotation=45)
	ax2 = Axis(f[2,2], xlabel = "energy [keV]", ylabel = "counts / 100keV" , title = "zoomed to 2500,3500", limits=(2500,3500, nothing,nothing), yaxisposition = :right, xticklabelrotation=45)

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
	- The peak of 0nu (in blue) is smeared quite a bit due to energy resolution + electron energy losses in gas. *In an ideal detector this should be delta function.*
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

# ╔═╡ 52e4668e-9804-492a-9deb-830f8cdedfcf
md"""
# 3. Model
"""

# ╔═╡ 68231ef5-9425-4e4c-98a8-586323938d45
md"""

"""

# ╔═╡ 1c44a17c-41a5-4ec2-a37e-b4c20e9407cc


# ╔═╡ 0033c567-5cf4-4067-bc50-00543dabc798


# ╔═╡ c69bf1d5-e7ae-4b31-8fe3-cbb872e68424
md"""
# To make the analysis fair when compring the two approaches I generate a sample spectrum from the presented background model:

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
	Since data is statistical process, it fluctuates a lot! 
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
	return data ./ 1000.0 
end;

# ╔═╡ 873c693a-23a9-4ec2-ac3a-04cd9c9ac40b
md"""
# From now on we will use this sample set to perform both analyses on.
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
\mathcal S(\bar{b}) \begin{cases} 
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
# From this map we can find the following values for the background model: 

- ### best ROI: ($( best_t12ESum[:minBinEdge] ) - $( best_t12ESum[:maxBinEdge] )) keV 
- signal efficiency in ROI: $(round(effbb, sigdigits=3) *100) %
- expected background counts $\bar{b}$ in ROI: $(round(expBkgESum, sigdigits=3) ) 
- corresponding $\mathcal S(\bar{b})$: $(round(nexcl, sigdigits=3) ) 
- ### max sensitivity: $(round(best_t12ESum[:maxBinCount], sigdigits=3)) yr

### Background overview table:
$(
pretty_isotopoes
)
"""

# ╔═╡ 6ea1a763-3268-4f77-83a0-f6065d77eb0a
md"""
## However, in order to be able to **compare** the sensitivity with Bayessian analysis, we have to calculate the sensitivity **NOT** from the background model but from the **sample data-set**. We do this in the following way: 

!!! note "Frequentist procedure for data-set"
	- First we cut the sample data-set to just energies close to the ROI (a reqion bit wider that ROI).
	- This data-set we use as a fitting region. I use an exponential fit here, from which we obtain $\lambda$ as the background's decay rate. 
	- We integrate the Exponential within the boundaries of the ROI to obtain $\bar{b} = (b-a) \int_{a}^{b} N e^{-\lambda E} dE$ (where a, b are the ROI boundaries; λ is decay constant and N is normalization).
	- We convert sample's $\bar{b}$ to $\mathcal S$.
	- Lastly, plug into sensitivity equation and calulate sample's sensitvitiy.
"""

# ╔═╡ 796d7a86-dae3-4e0c-9dc7-157d245425e9
md"""
a = $(@bind a PlutoUI.Slider(1.600:0.100:2.700, default=2.000; show_value=true)) MeV

b = $(@bind b PlutoUI.Slider(3.100:0.100:3.500, default=3.100; show_value=:true)) MeV
"""

# ╔═╡ 2485d73a-28bc-461a-be00-0083dba8a5fb
@bind clicked1 PlutoUI.Button("Generate data")

# ╔═╡ 1554ee0b-4a9c-49e4-badc-b000c767ed2d
data = let
	clicked
	clicked1
	generate_sample_data(h1d_full_bkg)
end

# ╔═╡ 23f06d08-e9ad-4de0-a152-fc5d47d6d103
h1d_sample = let 
	Hist1D(data; binedges=binningDict[:SumEMeV]) 
end;

# ╔═╡ 26f30acb-b8ab-4619-bdb2-b49107c0a317
let
	f = Figure(size = (1050, 600), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel = "energy [MeV]", ylabel = "counts / 100keV" , title = "sample background data \nlinear scale", xticklabelrotation=45)
	ax2 = Axis(f[1,2], xlabel = "energy [MeV]", ylabel = "counts / 100keV" , title = "log scale", yscale =log10, xticklabelrotation=45)
	
	st = hist!(ax, h1d_sample, label ="background")
	st2 = hist!(ax2, h1d_sample, label ="log background",)
	ylims!(ax2, 1e-1, 1e5)
	f
end

# ╔═╡ 3afc4c70-452c-4d7b-8795-cf82567d7475
function likelihood_a(n_B::Real, n_S::Real, Q::Real, sigma::Real, lambda::Real, x::Real)
    # res = (pdf(Exponential(lambda), x) + pdf(Normal(Q, sigma),  x)) 
	lambda <= 0.0 && return 1.0
    res = (n_B * pdf(Exponential(lambda), x) + n_S * pdf(Normal(Q, sigma),  x)) * inv( n_B + n_S ) 
	res = res <= 0.0 ? 1.0 : res
	return res
end;

# ╔═╡ be3221e1-553d-40d4-8ec1-c9c272be7251
function nll_exp(params, d::Vector{<:Real})
    n_B, n_S, Q, sigma, lambda = params
    ll = @. log(likelihood_a(n_B, n_S, Q, sigma, lambda, d))
    return -sum(ll)
end;

# ╔═╡ 0af4a737-778b-484d-a317-28c9f4b07668
BigFloat(1e20)

# ╔═╡ 5d6fb41a-8146-4423-94b2-dc2e78495a25
f_exp(x, p) = p[1]./exp.(x.*p[2])

# ╔═╡ 1208d6be-4cc6-4cb9-a16f-4e7323d38e59
initial_p = [float(1e10), 12.4]

# ╔═╡ a449ff98-c5cc-4f88-a700-edf2d610eab3
# ╠═╡ disabled = true
#=╠═╡
b_exp_sample = exp_integral(fit.param[2], fit.param[1], a,b) * 17.5 
  ╠═╡ =#

# ╔═╡ 095600d3-bf7d-4373-9e92-b92557010b49
sample_ROI_data = filter( x-> a <= x <= b, data )

# ╔═╡ c48ff5b5-cfcc-48b1-9859-70d5c0fae260


# ╔═╡ a8025e0c-a17a-467f-9873-6999847d184f
h1d_sample_ROI = restrict(h1d_sample, a, b)

# ╔═╡ 8eaffbf3-39f3-41c4-afbd-7619b9613062
bincenters(h1d_sample_ROI)

# ╔═╡ fb4c4cbf-1e6b-495c-a7d7-f75c0a15cce0
bincounts(h1d_sample_ROI)

# ╔═╡ 108822d9-e25d-4aaa-bb61-33f30faff50a
fit = curve_fit(f_exp, bincenters(h1d_sample_ROI), bincounts(h1d_sample_ROI), initial_p)

# ╔═╡ c5f99485-1581-4b3f-8d50-af2cabf8dd39
h1d_sample_ROI |> integral

# ╔═╡ e20a262d-362b-40d1-83a3-826525251f50
function chi_square(bincounts, f, xs; nparams= 1)
	ndf = length(bincounts) - nparams
	χ² = 0.0
	for i in eachindex(bincounts)
		if ( bincounts[i] >= 5 )
			χ² += (f(xs[i]) - bincounts[i])^2 / f(xs[i])
		else
			ndf -= 1
		end
	end
	
	return χ², ndf
end;

# ╔═╡ 8f81bf60-cb6a-4d24-a793-18cb621c1afb
let l = fit.param[2], N = fit.param[1], h = h1d_sample_ROI, fit_f(d) = f_exp(d, fit.param)
	
	f = Figure(size = (800, 500), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel ="energy", ylabel="counts", xticks= range(a, b, 5), yscale =log10)
	
	bin_edges = binedges(h)
	n = sum(bincounts(h))
	min_bin, max_bin = extrema(bin_edges) # used for scaling on the x-axis
	
	xs = range(min_bin, max_bin, length=100)  
	ys = fit_f.(xs) 

	chisq, ndf = chi_square( bincounts(h), fit_f, bincenters(h))

	p = plot!(ax, h, label = "data", xlabel ="energy", ylabel="normalized counts", )
	vspan!(ax, best_t12ESum[:minBinEdge]/1000.0, best_t12ESum[:maxBinEdge]/1000.0, color = (:orange, 0.3), label ="ROI" )
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

# ╔═╡ a849c9be-f4bf-45a6-8ba9-85224077532e
exp_integral(12.4, 5.8e14, best_t12ESum[:minBinEdge]/1000. , best_t12ESum[:maxBinEdge]/1000.) * SNparams["tYear"] * SNparams["foilMass"]* SNparams["a"]

# ╔═╡ 890ef7fc-2649-442e-beaf-8f9db8c69220
b_exp_sample = exp_integral(fit.param[2], fit.param[1], best_t12ESum[:minBinEdge]/1000. , best_t12ESum[:maxBinEdge]/1000.) * SNparams["tYear"] * SNparams["foilMass"]* SNparams["a"]

# ╔═╡ df989bf7-ee2c-4912-bee8-42afeeca70e1
T12_sample_freq = round(get_tHalf(SNparams, effbb, b_exp_sample, α), sigdigits = 4);

# ╔═╡ ef493c16-656c-4253-bb46-62e0aefae4c9
md"""
## From this fit we can extract expected number of background counts in ROI per 17.5kgy. 

### $\int_{ROI} \lambda e^{-\lambda E}dE = \bar{b} =$ $(round(b_exp_sample, sigdigits=4))

#### This value can be plugged into the sesnitivity equation to obtain: $T^{1/2} \geq$ $T12_sample_freq y

!!! danger "Note on the result:" 
	The value obtained here seems a bit **higher** than the one obtained from the model. The problem here is that with so few data-points in the region of interest, the numbers fluctuate a lot and thus the result is dominated by uncertainty. In the final calculation, once the acitivities are actually measured from data, the background model can be optimized and used more properly. This example is just a small showcasing. 
"""

# ╔═╡ 46086903-ec96-4f86-8081-ddfb713fde30
md"""
# 5. Bayessian approach: "Here's what the data tells me..."

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

	To build the model we look at the sample data:

"""

# ╔═╡ 28a5b0b3-8055-49a4-8623-688c785be07d
let h = (h1d_sample_ROI) 
	f = Figure(size = (800, 500), fontsize= 22, fonts = (; regular = "TeX"))
	ax = Axis(f[1,1], xlabel ="energy", ylabel="counts", xticks= range(a, b, 5))
	p = plot!(ax, h, label = "data", xlabel ="energy", ylabel="normalized counts", )
	f[1,2] = Legend(f, ax, merge = true )
	f
end


# ╔═╡ 2c04e398-8c31-4f9b-a00a-528c9bf60358
md"""

!!! update "Likelihood, prior, posterior...:"
	As with any attempts at making a model with the hopes of describing data, we must make a number of assumption:
	1. Data consists of some amount of background $n_b$ and some amount of signal $n_s$ (i.e. the Signal + Background model)
	2. Background is distributed exponentially with decay constant $\lambda$: $pdf(b) = \lambda e^{-\lambda E}$ 

"""

# ╔═╡ 0dc66ff5-c2c5-4b86-a971-f9f01e7dc7ff
md"""
# The end! 😏
"""

# ╔═╡ aa75ac59-2e73-440e-90ed-ff397c696892
@bind mu PlutoUI.Slider(5:15; show_value=true)

# ╔═╡ 5a7f967c-32f4-448e-8382-39547104a499
f, ax, p = plot(Normal(mu, 20),);

# ╔═╡ 62a82632-ef0d-4617-8b3d-981ab675d56a
xlims!(ax, -100, 100)

# ╔═╡ 7d7e0703-12d1-4f37-8998-f966a7af6e62
f

# ╔═╡ Cell order:
# ╠═e4c0246e-b67a-4776-8f15-6efd5aec1f00
# ╠═3bd341f6-eae9-41f5-b5eb-b6c3e76cbaf2
# ╠═e0d5731a-086d-11ef-0315-5705b559453d
# ╠═2c1253c6-5601-4701-944c-03fcad58c175
# ╠═c56b1bac-c440-4f71-a7a3-088b3042dfad
# ╠═099b4434-ca9d-4748-ba3e-7f097e3229e3
# ╠═773e7c01-dc0f-464a-9696-9a89e2302aa8
# ╠═4b9a7ce2-3e5f-4995-af22-985b5fb0f14e
# ╠═acdec082-7a61-489f-9642-9a459ab9f04a
# ╠═503673d7-7dac-43eb-a493-0aea57674720
# ╠═50e7af06-2bce-4914-a867-ec6261a6714b
# ╠═bee405ae-e2bd-48b8-8de8-652b99bbd8c4
# ╠═955810a1-9ef6-4531-8156-02b070c29676
# ╠═11b39aa3-2c3b-4d5f-8ef0-e938c3cb3399
# ╠═e4c6c2b0-2afe-4667-894c-71c8476239c6
# ╠═94ad3e49-34d8-45e5-ac65-a9241136f226
# ╠═875130d7-f0e3-4be0-a963-c427712ccdcf
# ╠═557e5713-6563-4e4e-b0ce-a9ee67cf12e6
# ╠═b1991653-1363-4810-a78d-923c4ce64daa
# ╠═a65e5feb-fe4a-41a0-9cc9-530c0c90a317
# ╠═21b41886-2a18-40d2-b57c-62862b30b8fe
# ╠═cede3473-fae8-4998-ba6e-01f9d05ad54c
# ╠═6821c511-ef71-401e-9f36-1719886e7b54
# ╠═52165d49-4a46-41d5-9803-20047b1ce12e
# ╠═3cbcfaf5-b367-40e2-9093-d10f852d403a
# ╠═f57a6ed8-2332-4ab1-abef-27c44881c7af
# ╠═84c47066-c7b7-4653-8e7a-9d24305fe091
# ╠═ec440471-bd5a-44ca-bd8c-7be723c07179
# ╠═c06c24af-9387-4094-b8f5-ac89e01230c0
# ╠═ce3b0244-4acf-4670-bfca-9d9715ee0d93
# ╠═8d9a138e-e015-411c-bc27-966d06fdb901
# ╠═3c5fadf4-5791-4f33-9ce1-6394611b977b
# ╠═54824dab-3172-46f2-b130-31f55a8262fb
# ╠═bceec9c1-785e-4aec-8195-82d1db4ee838
# ╠═f5f2503a-1809-49c2-8839-3ea17c99b286
# ╠═84acc905-d854-42c2-a862-e362410205b2
# ╠═52e4668e-9804-492a-9deb-830f8cdedfcf
# ╠═68231ef5-9425-4e4c-98a8-586323938d45
# ╠═1c44a17c-41a5-4ec2-a37e-b4c20e9407cc
# ╠═0033c567-5cf4-4067-bc50-00543dabc798
# ╠═c69bf1d5-e7ae-4b31-8fe3-cbb872e68424
# ╠═58241414-5ec2-463a-8422-f7bd8e30b606
# ╠═91546e37-2d22-4616-9d19-a4d973e80e54
# ╠═d5b67056-4212-40b6-b25c-007cacd0509a
# ╠═26f30acb-b8ab-4619-bdb2-b49107c0a317
# ╠═06e1024a-5afb-4cd8-a75e-4bcfb4a5ab8f
# ╠═1554ee0b-4a9c-49e4-badc-b000c767ed2d
# ╠═23f06d08-e9ad-4de0-a152-fc5d47d6d103
# ╠═a47c9336-33de-4c89-bc64-3e1e12036243
# ╠═873c693a-23a9-4ec2-ac3a-04cd9c9ac40b
# ╠═78c184cf-e6c2-485d-9876-e3811ae00c8f
# ╠═6ed9cdf8-6d18-456e-934f-035f8e2df1d3
# ╠═883a548d-54ec-4242-9e65-38118045a517
# ╠═9be29518-9b98-4f6d-a504-40a27ada1107
# ╠═e5cba772-d3a1-4b5b-9897-164183366e04
# ╠═7e9ac44f-b1af-477b-a9ab-762e12a9d294
# ╠═bb31751d-eaa8-4b7e-aad0-cf894a6ee7f7
# ╠═14d69080-2f1b-48a6-8c05-d58b5b4d3992
# ╠═196c8517-cb7a-4286-b630-98488d8bbe18
# ╠═abc26356-b76e-4bc0-a1fd-c9aed18c0841
# ╠═be30be54-ccad-4e7e-a758-f55b4dfbfa35
# ╠═a4b6a249-d490-4c97-af08-aaa7548d0b86
# ╠═6879f6a2-6d14-45de-99d6-f0bf0baceb85
# ╠═a859ecc0-bfae-4c5b-87da-5d5804a60ba1
# ╠═010d6d22-df4c-4dbf-b31e-21629201b8c6
# ╠═cd9e7bcf-193c-4c2d-9e42-5fd9f0004886
# ╠═ff1cc95d-bb3f-49b0-b988-8cdeca763b9d
# ╠═6356676f-32a7-452d-b6ee-3d6267c9e9fa
# ╟─6ea1a763-3268-4f77-83a0-f6065d77eb0a
# ╟─796d7a86-dae3-4e0c-9dc7-157d245425e9
# ╠═2485d73a-28bc-461a-be00-0083dba8a5fb
# ╟─8f81bf60-cb6a-4d24-a793-18cb621c1afb
# ╟─ef493c16-656c-4253-bb46-62e0aefae4c9
# ╠═3afc4c70-452c-4d7b-8795-cf82567d7475
# ╠═be3221e1-553d-40d4-8ec1-c9c272be7251
# ╠═039334ac-d249-4c5f-bc97-e1180b52b47b
# ╠═0af4a737-778b-484d-a317-28c9f4b07668
# ╠═5d6fb41a-8146-4423-94b2-dc2e78495a25
# ╠═1208d6be-4cc6-4cb9-a16f-4e7323d38e59
# ╠═8eaffbf3-39f3-41c4-afbd-7619b9613062
# ╠═fb4c4cbf-1e6b-495c-a7d7-f75c0a15cce0
# ╠═108822d9-e25d-4aaa-bb61-33f30faff50a
# ╠═df989bf7-ee2c-4912-bee8-42afeeca70e1
# ╠═a449ff98-c5cc-4f88-a700-edf2d610eab3
# ╠═a849c9be-f4bf-45a6-8ba9-85224077532e
# ╠═890ef7fc-2649-442e-beaf-8f9db8c69220
# ╠═c5f99485-1581-4b3f-8d50-af2cabf8dd39
# ╠═095600d3-bf7d-4373-9e92-b92557010b49
# ╠═c48ff5b5-cfcc-48b1-9859-70d5c0fae260
# ╠═a8025e0c-a17a-467f-9873-6999847d184f
# ╠═e20a262d-362b-40d1-83a3-826525251f50
# ╠═79b85562-16a8-423d-ac7d-16d9ab89c61b
# ╠═46086903-ec96-4f86-8081-ddfb713fde30
# ╠═05c42c02-7306-412b-b205-74701de21634
# ╠═e54eb759-9274-4df3-9cc3-12ba8925b067
# ╠═339077ec-a94d-4d36-9003-3dc690cdc00a
# ╠═d941baf4-d60a-40f8-a35d-bf76e16dbeda
# ╠═742f7fd5-01e7-4731-a7a1-ad84bb86b8c7
# ╠═278988c0-6da6-46f2-a50d-89ea32edb69a
# ╠═4cd71a16-3bb8-4d81-b1f7-f67eb6cad196
# ╠═28a5b0b3-8055-49a4-8623-688c785be07d
# ╠═2c04e398-8c31-4f9b-a00a-528c9bf60358
# ╠═0dc66ff5-c2c5-4b86-a971-f9f01e7dc7ff
# ╠═a5b84247-8464-4a80-af07-5309d7bcaef9
# ╠═5a7f967c-32f4-448e-8382-39547104a499
# ╠═62a82632-ef0d-4617-8b3d-981ab675d56a
# ╠═aa75ac59-2e73-440e-90ed-ff397c696892
# ╠═7d7e0703-12d1-4f37-8998-f966a7af6e62
