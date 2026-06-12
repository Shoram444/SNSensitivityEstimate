using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
using DataFramesMeta, CSV
using CairoMakie
using StatsPlots


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
h_tot = h_signal + h_bkg 

StatsPlots.plot(normalize(h_signal, width = false),label = "Signal", color = :dodgerblue3, fillalpha = 0.3)
StatsPlots.plot!(normalize(h_bkg, width = false),  label = "Background", color = :orange, fillalpha = 0.3)
StatsPlots.plot!(normalize(h_tot, width = false), label = "Total", color = :green, fillalpha = 0.3)


# Generate sample spectra 
pseudo_bkg = StatsBase.sample(binedges(h_bkg)[1:end-1], Weights(bincounts(h_bkg)), nb_pass)
pseudo_signal = StatsBase.sample(binedges(h_signal)[1:end-1], Weights(bincounts(h_signal)), 10)
data_hist = Hist1D(; binedges = 0:0.1:10)

push!.(data_hist, pseudo_bkg)
push!.(data_hist, pseudo_signal)

StatsPlots.plot(normalize(data_hist, width = false), label = "Pseudo-data", color = :black, fillalpha = 0.3)

bkg_hist_normed = [normalize(h_bkg, width = true)]
signal_hist_normed = normalize(h_signal, width = true)

my_likelihood = make_hist_likelihood_uniform(
    data_hist,
    signal_hist_normed,
    bkg_hist_normed
)


# Define priors 
prior = NamedTupleDist(
    As = Uniform(0, 2),  
    Ab = [Uniform(0, 2)],
)   

# Define sampling options
burnin = MCMCMultiCycleBurnin(max_ncycles = 50, nsteps_final=1000)
mcmcalgo = RandomWalk()

posterior = PosteriorMeasure(my_likelihood, prior)
samples, _ = bat_sample( posterior, MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 5*10^4, nchains = 4))

# plot posterior
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


# --- BAT fit plotting helpers -------------------------------------------------

function _extract_As_Ab(binned_unshaped_samples)
    As = Float64[]
    Ab = Float64[]

    for par in binned_unshaped_samples.v
        push!(As, par[1])

        ab_raw = par[2]
        if ab_raw isa Number
            push!(Ab, float(ab_raw))
        else
            push!(Ab, float(first(ab_raw)))
        end
    end

    return As, Ab
end

function _model_curves(binned_unshaped_samples, data_hist, signal_hist_normed, bkg_hist_normed; n_grid = 300)
    As, Ab = _extract_As_Ab(binned_unshaped_samples)

    bw = inv(binedges(data_hist).inv_step)
    xlo = first(binedges(data_hist)) + 0.5 * bw
    xhi = last(binedges(data_hist)) - 0.5 * bw
    x = range(xlo, stop = xhi, length = n_grid)
    ndata = integral(data_hist)
    n_samples = length(As)

    ymat = Matrix{Float64}(undef, n_samples, n_grid)
    bkg_template = bkg_hist_normed isa AbstractVector ? first(bkg_hist_normed) : bkg_hist_normed

    for i in 1:n_samples
        h = As[i] * signal_hist_normed + Ab[i] * bkg_template
        y = map(xi -> begin
            yi = lookup(h, xi)
            ismissing(yi) ? 0.0 : float(yi)
        end, x)
        @views ymat[i, :] .= y .* ndata .* bw
    end

    return x, ymat
end

function plot_fit_samples(
    binned_unshaped_samples,
    data_hist,
    signal_hist_normed,
    bkg_hist_normed;
    n_plot_samples = 100,
    n_grid = 300,
)
    x, ymat = _model_curves(
        binned_unshaped_samples,
        data_hist,
        signal_hist_normed,
        bkg_hist_normed;
        n_grid,
    )

    n_total = size(ymat, 1)
    n_take = min(n_plot_samples, n_total)
    idx = round.(Int, range(1, n_total, length = n_take))

    bin_mids = midpoints(binedges(data_hist))
    data_counts = bincounts(data_hist)
    data_err = sqrt.(max.(data_counts, 1.0))

    fig = Figure(size = (860, 560))
    ax = Axis(
        fig[1, 1],
        xlabel = "Energy [a.u.]",
        ylabel = "Counts / bin",
        title = "BAT fit: posterior sample overlays",
        xgridvisible = false,
        ygridvisible = true,
    )

    band!(ax, [first(x), last(x)], [0.0, 0.0], [maximum(data_counts) * 1.15, maximum(data_counts) * 1.15],
        color = (:gray95, 0.35))

    for i in idx
        lines!(ax, x, @view(ymat[i, :]), color = (:dodgerblue3, 0.14), linewidth = 1.2)
    end

    errorbars!(ax, bin_mids, data_counts, data_err, color = :black, whiskerwidth = 8)
    CairoMakie.scatter!(ax, bin_mids, data_counts, color = :black, markersize = 6)

    CairoMakie.xlims!(ax, first(x), last(x))
    fig
end

function plot_fit_best_with_bands(
    binned_unshaped_samples,
    data_hist,
    signal_hist_normed,
    bkg_hist_normed;
    n_grid = 300,
)
    x, ymat = _model_curves(
        binned_unshaped_samples,
        data_hist,
        signal_hist_normed,
        bkg_hist_normed;
        n_grid,
    )

    q025 = vec(mapslices(v -> quantile(v, 0.025), ymat, dims = 1))
    q16 = vec(mapslices(v -> quantile(v, 0.16), ymat, dims = 1))
    q50 = vec(mapslices(v -> quantile(v, 0.50), ymat, dims = 1))
    q84 = vec(mapslices(v -> quantile(v, 0.84), ymat, dims = 1))
    q975 = vec(mapslices(v -> quantile(v, 0.975), ymat, dims = 1))

    bin_mids = midpoints(binedges(data_hist))
    data_counts = bincounts(data_hist)
    data_err = sqrt.(max.(data_counts, 1.0))

    fig = Figure(size = (860, 560))
    ax = Axis(
        fig[1, 1],
        xlabel = "Energy [a.u.]",
        ylabel = "Counts / bin",
        title = "BAT fit: best estimate with 68% and 95% bands",
        xgridvisible = false,
        ygridvisible = true,
    )

    band!(ax, x, q025, q975, color = (:lightskyblue, 0.30), label = "95% band")
    band!(ax, x, q16, q84, color = (:deepskyblue4, 0.32), label = "68% band")
    lines!(ax, x, q50, color = :navy, linewidth = 3, label = "Median fit")

    errorbars!(ax, bin_mids, data_counts, data_err, color = :black, whiskerwidth = 8)
    CairoMakie.scatter!(ax, bin_mids, data_counts, color = :black, markersize = 6, label = "Data")

    axislegend(ax, position = :rt)
    CairoMakie.xlims!(ax, first(x), last(x))
    fig
end

fig_samples = plot_fit_samples(
    binned_unshaped_samples,
    data_hist,
    signal_hist_normed,
    bkg_hist_normed;
    n_plot_samples = 100,
)

fig_bands = plot_fit_best_with_bands(
    binned_unshaped_samples,
    data_hist,
    signal_hist_normed,
    bkg_hist_normed,
)
