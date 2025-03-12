import Pkg
Pkg.activate(".")

using Random, LinearAlgebra, Statistics, Distributions, StatsBase, Plots, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist
using SpecialFunctions

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

plot(samples, :(a), nbins=50)
