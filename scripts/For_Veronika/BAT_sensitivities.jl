using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
using DataFramesMeta, CSV
using CairoMakie

# first we generate mock data
# we assume gaussian signal process and exponential background
# we assume a dataset of 100 signal events, nS and 1000 background nB 
# since the estimate relies on the efficiency of signal detection 
# (which is calculated as the number of events that pass analysis counts
# np over the total number of simulated events ns; ε = np/ns, we must
# make up some number for ns > nS.

ns_sim = 10_000_000 # mock simulated number of signal events (before datacuts) 
nb_sim = 10_000_000 # mock simulated numbed of background events (before datacuts)


ns_pass = 1_000_000 # n signal events (after data-cuts)
nb_pass = 10_000 # n backgoround events (after data-cuts)

es = ns_pass/ns_sim # signal efficiency
eb = nb_pass/nb_sim # background efficiency

mu_signal = 3 # signal mean value - Q-value keV
sigma_signal = 0.3*mu_signal # signal sigma - energy resolution keV

lambda_background = 3 # decay constant of background

n_signal_dist = rand(Normal(mu_signal, sigma_signal), ns_pass)
n_background_dist = rand(Exponential(lambda_background), nb_pass)

# plot the distributions
h_signal = Hist1D(; binedges = 0:0.1:10)
h_bkg = Hist1D(; binedges = 0:0.1:10)

push!.(h_signal, n_signal_dist)
push!.(h_bkg, n_background_dist)
h_tot = (h_signal + h_bkg )|> normalize

stephist(h_signal)
stephist!(h_bkg)
stephist!(h_tot)
current_figure()

# Generate sample spectra 
pseudo_bkg = StatsBase.sample(binedges(h_bkg)[1:end-1], Weights(bincounts(h_bkg)), nb_pass)
pseudo_signal = StatsBase.sample(binedges(h_signal)[1:end-1], Weights(bincounts(h_signal)), 10)
data_hist = Hist1D(; binedges = 0:0.1:10)

push!.(data_hist, pseudo_bkg)
push!.(data_hist, pseudo_signal)

bkg_hist_normed = [normalize(h_bkg, width = true)]
signal_hist_normed = normalize(h_signal, width = true)

my_likelihood = make_hist_likelihood_uniform(
    data_hist,
    signal_hist_normed,
    bkg_hist_normed
)


# Define priors 
prior = NamedTupleDist(
    As = Uniform(0, 1),  
    Ab = [Uniform(0, 1)],
)   

# Define sampling options
burnin = MCMCMultiCycleBurnin(max_ncycles = 50, nsteps_final=1000)
mcmcalgo = RandomWalk()

posterior = PosteriorMeasure(my_likelihood, prior)
samples, _ = bat_sample( posterior, MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 10^4, nchains = 4))

# plot posterior
using StatsPlots
StatsPlots.plot(samples)

samples_df = DataFrame(samples)


# Fit to sensitivity
binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
nDataPoints = integral(data_hist)
muS = [par[1] for par in binned_unshaped_samples.v]

@show exp_mu_signal_90 = BAT.smallest_credible_intervals(muS, nsigma_equivalent= 1.65)[1].right * nDataPoints

Na = SNparams["Nₐ"]
m = SNparams["foilMass"] * SNparams["a"]
t = SNparams["tYear"]
W = SNparams["W"]

es 

Thalf = log(2) * (Na * m * t * es / W) / exp_mu_signal_90


Plots.plot(samples, :As)

