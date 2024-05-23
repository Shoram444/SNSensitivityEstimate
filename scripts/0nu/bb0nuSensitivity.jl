using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Distributions
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

all_processes = load_processes("fal5_8perc", "sumE")

signal = get_process("bb0nu_foil_bulk", all_processes)
background = [
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

Q_keV = SNparams["Q"]
Q_MeV = Q_keV / 1000.0
α = 1.64485362695147


t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="formula")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)
effbb = lookup(signal, best_t12ESum)
ThalfbbESum = round(get_tHalf(SNparams, effbb, expBkgESum, α), sigdigits=3)

lbl = "$(best_t12ESum[:minBinEdge]) - $(best_t12ESum[:maxBinEdge]) keV 
      b  = $(round(expBkgESum, sigdigits = 3)) 
      T12 ≥  $(round(ThalfbbESum, sigdigits=2)) yr 
      ε = $(round(effbb, sigdigits = 3)*100)%"

let 
    f = Figure(size=(600, 400))
    a = Axis(f[1,1], xlabel = "min_ROI [keV]", ylabel = "max_ROI [keV]")
    p = plot!(a, t12MapESum)
    text!(a, 2000, 500, text=lbl)
    Colorbar(f[1,2], p, label="sensitivity [yr]", scale=log10)
    f
end

t = range(0, 6, 60)

t12(t, e, b) = get_tHalf(
    SNparams["W"],
    SNparams["foilMass"],
    SNparams["Nₐ"],
    t,
    SNparams["a"],
    e,
    b / SNparams["tYear"] * t,
    α;
    approximate="formula"
)

let 
    f = Figure(size=(600, 400))
    a = Axis(f[1,1], xlabel = "detector life-time [yr]", ylabel = "sensitivity [yr]", limits= (0,6, nothing, nothing))
    p = lines!(a, t, t12.(t, effbb,expBkgESum), label = "FC: 8% resolution; simulation")
    axislegend(a, position = :lt)
    f
end

a = 2500.0
b = 3500.0

ROI_a, ROI_b = best_t12ESum[:minBinEdge], best_t12ESum[:maxBinEdge] 

h1d_background = get_bkg_counts_1D.(background) |> sum
h1d_background = restrict(h1d_background, a, b)

pseudo_data = generate_pseudo_data(h1d_background)  


### Fitting
using Optim, Turing, MCMCChains
model1 = exponential_normal_mixture(pseudo_data .- a, a, b)

mle_fit = optimize(model1, MLE(), NelderMead()).values
chains = Turing.sample(model1, NUTS(), 10_000)

function fit_function_exponential_normal_model( lambda, sigma, nSig, nBkg, x, Q, a, b )
    nSig * pdf(Normal(Q, sigma), x-a ) + nBkg*pdf( Exponential(lambda), x-a )
end

xs = range(a,b,100)
sigma_at_Q = get_sigma_keV(Q_keV, 0.08)
# mle fit params
mu_mle, sigma_mle, l_mle, ns_mle, nb_mle = mle_fit
ns_mle = ns_mle * length(pseudo_data)
nb_mle = nb_mle * length(pseudo_data)
# Bayess fit params
summary_chains = summarize(chains) |> DataFrame
mu_bayes, sigma_bayes,l_bayes, ns_bayes, nb_bayes = summary_chains[:, 2]
ns_bayes = ns_bayes * length(pseudo_data)
nb_bayes = nb_bayes * length(pseudo_data)


ys_bayess = fit_function_exponential_normal_model.(l_bayes, sigma_at_Q, ns_bayes, nb_bayes, xs, Q_keV, a,b) .* step(binningDict[:SumE]) 
ys_mle = fit_function_exponential_normal_model.(l_mle, sigma_at_Q, ns_mle, nb_mle, xs, Q_keV, a,b) .* step(binningDict[:SumE]) 

let 
    f = Figure(size=(600, 400))
    a = Axis(f[1,1], xlabel = "energy [keV]", ylabel = "counts")
    p = plot!(a, h1d_background, label ="pseudo-data")
    lines!(a, xs, ys_mle, label = "fit_mle", color=Makie.wong_colors()[2], linewidth=2, linestyle=:dash)
    lines!(a, xs, ys_bayess, label = "fit_bayes", color=Makie.wong_colors()[3], linewidth=2, linestyle=:dot)
    axislegend(a, position = :rt)
    f
end

credible_interval_90 = hpd(chains; alpha=0.1) |> DataFrame
mu_U = credible_interval_90[4,3] * length(pseudo_data)

using QuadGK
f(x) = fit_function_exponential_normal_model(l_mle, sigma_at_Q, 0, nb_mle, x, Q_keV, a, b) 
b_exp_mle, error = quadgk(f, ROI_a, ROI_b)

function get_tHalf_bayess( SNparams, eff, mu_U )
	return log(2)*SNparams["Nₐ"]*SNparams["foilMass"]*SNparams["a"]*SNparams["tYear"]/SNparams["W"]*eff/mu_U 
end

let 
    f = Figure(size=(650, 450), fontsize = 16)
    a = Axis(
        f[1,1], 
        xlabel = "detector life-time [yr]", 
        ylabel = "sensitivity [yr]", 
        limits= (0,6, 0, nothing), 
        title= "SuperNEMO sensitivity as a function of detector life-time \ncalculated at 90% confidence level ")
    p1 = lines!(a, t, t12.(t, effbb,expBkgESum), label = "simulation", linewidth=3)
    p2 = lines!(a, t, t12.(t, effbb,b_exp_mle), label = "pseudo-data", linewidth=3)
    s1 = scatter!(a, 2.88, get_tHalf_bayess(SNparams, effbb,2.37),  color = :green, marker=:star6, markersize= 16)
    h1 = hlines!(a, [4.6e24], color = :black,  linewidth=3)
    vlines!(a, [2.88], color=(:grey, 0.5), linestyle=:dot, linewidth=2)
    text!(
        3, 3.2e24,
        text = "Preliminary",
        rotation = pi/6,
        align = (:center, :baseline),
        color = (:red, 0.3),
        fontsize = 40
    )
    axislegend(a,[p1,p2], ["simulation", "pseudo-data"], "frequentist 90% CL", position = :lt, margin= (15,0,0,15))
    axislegend(a, [s1], ["pseudo-data"],"Bayes 90% CI\n17.5 kg.y exposure", position = :lt, margin= (15,0,0,115), padding=(7,7,3,3), tellwidth = false,)
    axislegend(a, [h1], ["Cupid-0 data"],"Bayes 90% CI\n8.82 kg.y exposure", position = :rb, margin= (0,15,15,0), padding=(7,7,3,3), tellwidth = false,)
    save("SuperNEMO_sensitivity.png", f)
    f
end



using BAT, DensityInterface

prior = distprod(
    nB=Uniform(1e-5, 1e3), #Uniform(1e-5, 45),
    nS=Uniform(1e-15, 10), #Uniform(1e-5, 15),
    lambda = Uniform(1e-5,1e3)
)

unbinned_posterior = PosteriorMeasure(BAT_normal_exponential_model(pseudo_data, Q_keV, 0.08), prior)
unbinned_samples = bat_sample(unbinned_posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result
unbinned_unshaped_samples, f_flatten = bat_transform(Vector, unbinned_samples)
muSsUnbinned = [m[2] for m in unbinned_unshaped_samples.v]

MCMCChains._hpd(muSsUnbinned; alpha=0.1)

