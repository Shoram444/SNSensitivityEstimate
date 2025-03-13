using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())

using SensitivityModule
using Random, LinearAlgebra, Statistics, Distributions, Plots, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions

mu1 = 1000
lambda = 0.001
sigma1 = sqrt(mu1)
n1 = 1
n2 = 1000
n = n1 + n2
ε = 0.1
d1 = rand(Normal( mu1, sigma1), n1) .+ ε
d2 = rand(Uniform(0, 3000), n2) .+ ε

data = vcat(
    d1, d2
)


h1 = normalize(Hist1D(d1; binedges = 0:100:3000), width = true)
h2 = normalize(Hist1D(d2; binedges = 0:100:3000), width = true)



data_hist = Hist1D(data; binedges = 0:100:3000)

true_par_values = (a = n1/sum(bincounts(data_hist)), ε = 0.1)

plot(data_hist)

my_pdf(h::Hist1D, x::Real) = ismissing(lookup(h, x)) ? 0.0 : lookup(h, x)

function f(par::NamedTuple{(:a, :ε)}, x::Real)
    total_rate = sum(par.a)
    a1 = par.a[1]
    a2 = 1.0 - a1
    th = normalize(a1*h1 + a2*h2 , width = true)
    return my_pdf(th, x) 
end


prior = distprod(
    a = Uniform(0, 1),
    ε = 0.1
)

function log_pdf_poisson(λ::T, k::U) where {T<:Real,U<:Real}
    R = float(promote_type(T,U))
    if λ >= 0 && k >= 0 && isinteger(k)
        a = iszero(k) ? R(k) : R(log(λ) * k)
        lg = logabsgamma(k + 1)[1]
        R(a - λ - lg)
    else
        R(-Inf)
    end
end

function make_hist_likelihood(h::Hist1D, f::Function)
    logfuncdensity(function(p)
        # Histogram counts for each bin as an array:
        counts = bincounts(h)
        n = sum(counts)
        # Histogram binning, has length (length(counts) + 1):
        binning = binedges(h)

        # sum log-likelihood value over bins:
        ll_value::Float64 = eps(Float64)
        
        for i in eachindex(counts)

            # Get information about current bin:
            bin_left, bin_right = binning[i], binning[i+1]
            bin_width = bin_right - bin_left
            bin_center = (bin_right + bin_left) / 2

            observed_counts = counts[i]

            # Simple mid-point rule integration of fit function `f` over bin:
            expected_counts = bin_width * f(p, bin_center) * n
            expected_counts == 0.0 && continue

            log_pdf_poisson(expected_counts, observed_counts) == -Inf && continue
            # Add log of Poisson probability for bin:
            ll_value += log_pdf_poisson(expected_counts, observed_counts)
            
        end
        ll_value
        return ll_value
    end )
end

my_likelihood = make_hist_likelihood(data_hist, f)

posterior = PosteriorMeasure(my_likelihood, prior)
samples, evals = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4))

marginal_modes = bat_marginalmode(samples).result

bat_report(samples)
true_par_values
plot(samples, :(a), nbins=50)

###################################
include(scriptsdir("Params.jl"))

# Load all the processes in the directory. Function `load_processes` takes two arguments:
# 1. dir::String -> the name of the directory where the root files are stored
# 2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
all_processes = load_data_processes("fal5_8perc_Boff_TIT_twoDistinct_edep_bcu", "sumE")

# declare which process is signal
signal = get_process("bb0nu_foil_bulk", all_processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
    get_process("K40_foil_bulk", all_processes),
    get_process("Pa234m_foil_bulk", all_processes),
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_nTotalSim!( signal, 0.98e8 )
# set_nTotalSim!( signal, 1e8 )

set_signal!(background[1], false)

# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 0.99e8 )
set_nTotalSim!( background[2], 0.96e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 0.76e8 )
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )


println("Processes initialized.")

#
Q_keV = SNparams["Q"]
α = 1.64485362695147


t12MapESum = get_tHalf_map(SNparams, α, signal, background...; approximate ="table")
best_t12ESum = get_max_bin(t12MapESum)
expBkgESum = get_bkg_counts_ROI(best_t12ESum, background...)

function get_sens_bayes(background, signal)
    ROI_a, ROI_b = 0, 3400

    bkg_hist = [(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)] |> sum
    bkg_hist_normed = normalize(bkg_hist, width = true)
    signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b), width = true)

    data_bkg = [first(FHist.sample(bkg_hist)) for i=1:rand(Poisson(round(Int, integral(bkg_hist))))] 
    data_hist = Hist1D( data_bkg; binedges= binedges(bkg_hist_normed) )

    function f1(par::NamedTuple{(:a, :ε)}, x::Real)
        total_rate = sum(par.a)
        a1 = par.a[1]
        a2 = 1.0 - a1
        th = normalize(a1*signal_hist_normed + a2*bkg_hist_normed , width = true)
        return my_pdf(th, x) 
    end

    my_likelihood = make_hist_likelihood(data_hist, f1)

    prior = distprod(
        a = Uniform(1e-20, 1),
        ε = 0.1
    )

    posterior = PosteriorMeasure(my_likelihood, prior)
    samples, evals = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^4, nchains = 4))

    marginal_modes = bat_marginalmode(samples).result

    binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
    nDataPoints = integral(data_hist)
    muS = [a[1] * nDataPoints for a in binned_unshaped_samples.v]

    """
        Stupid integration for getting upper limit - requires monotonely decreasing samples in vector form!
    """
    function get_interval_upper(data, CL; nbins=100)
        h = fit(Histogram, data; nbins = nbins)
        cs = cumsum(h.weights) ./ sum(h.weights)

        uppID = findfirst(x -> x >= CL, cs)
        midpoints(h.edges[1])[Int(uppID)]
    end

    exp_mu_signal_90 = get_interval_upper(muS, 0.9)
    Na = 6.02214e23
    m = 6.067
    t = 2.88
    W = 0.08192
    eff= lookup(signal, ROI_a, ROI_b)
    Thalf = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90
end

t = Float64[]
for i=1:100
    println(i)
    push!(t, get_sens_bayes(background, signal))
end

plot(t, st=:histogram, nbins = 10, xlabel = "Bayes sensitivity (yr)", label = "sample sensitivity")
vline!([median(t)], label = "Median = $(round(median(t), sigdigits = 3)) yr", color = :red, linewidth = 4)
vline!( [mean(t)], label = "Mean = $(round(mean(t), sigdigits = 3)) yr", color = :green, linewidth = 4)

plot(samples)
# plot(signal_hist_normed)

median(t)

scatter(t)




ROI_a, ROI_b = 0, 3400

bkg_hist = [(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)] |> sum
bkg_hist_normed = normalize(bkg_hist, width = true)
signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b), width = true)

data_bkg = [first(FHist.sample(bkg_hist)) for i=1:rand(Poisson(round(Int, integral(bkg_hist))))] 
data_hist = Hist1D( data_bkg; binedges= binedges(bkg_hist_normed) )

function f1(par::NamedTuple{(:a, :ε)}, x::Real)
    total_rate = sum(par.a)
    a1 = par.a[1]
    a2 = 1.0 - a1
    th = normalize(a1*signal_hist_normed + a2*bkg_hist_normed , width = true)
    return my_pdf(th, x) 
end

my_likelihood = make_hist_likelihood(data_hist, f1)

prior = distprod(
    a = Uniform(1e-20, 1),
    ε = 0.1
)

posterior = PosteriorMeasure(my_likelihood, prior)
samples, evals = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^4, nchains = 4))

marginal_modes = bat_marginalmode(samples).result

binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
nDataPoints = integral(data_hist)
muS = [a[1] * nDataPoints for a in binned_unshaped_samples.v]

"""
    Stupid integration for getting upper limit - requires monotonely decreasing samples in vector form!
"""
function get_interval_upper(data, CL; nbins=100)
    h = fit(Histogram, data; nbins = nbins)
    cs = cumsum(h.weights) ./ sum(h.weights)

    uppID = findfirst(x -> x >= CL, cs)
    midpoints(h.edges[1])[Int(uppID)]
end

exp_mu_signal_90 = get_interval_upper(muS, 0.9)

xs = 50:100:3350
ys = [f1((a = exp_mu_signal_90 / integral(data_hist), ε = 0.1), x) for x in xs]
plot(data_hist, label = "data",fillrange = 1e-5, fillcolor = :match)
plot!(xs, ys .* integral(data_hist) .* 100, linewidth = 4, ylims = (0, 1e4), label = "fit")
# plot!(xlims = (2000, 3200), ylims = (1e-5, 100))

ys_signal = bincounts(signal_hist_normed) .* 100 .* exp_mu_signal_90
plot!(xs, ys_signal, linewidth = 4, label = "signal")
plot!(yscale = :log10, ylims = (1.7e-1, 1e4))
savefig("scripts/0nu/Bayes_hist_models/fit.svg")