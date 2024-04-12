using BAT, Distributions, StatsPlots, DataFrames, StatsBase, LinearAlgebra
using DensityInterface, IntervalSets

bkgCounts_true = 20
sigCounts_true = 0
nTot = bkgCounts_true + sigCounts_true
Q = 2990
sigmaTrue = 130
minE = 2600
maxE = 3200
deltaE = 100

data = vcat(
    rand(Uniform(minE, maxE), bkgCounts_true),
    rand(Normal(Q, sigmaTrue), sigCounts_true)
)


function likelihood_a(p, x; minE=minE, maxE=maxE, Q=Q)
    (p.muB * pdf(Uniform(minE, maxE), x)  +
     p.muS * pdf(Normal(Q, p.sigma),  x)) *
    inv( p.muB + p.muS ) 
end

function likelihood_b(p, n)
    pdf(Poisson(p.muB + p.muS), n)
end

"""
    Stupid integration for getting upper limit - requires monotonely decreasing samples in vector form!
"""
function get_interval_upper(data, CL; nbins=100)
    h = fit(Histogram, data, nbins)
    cs = cumsum(h.weights)

    uppID = findfirst(x -> x >= CL, cs)
    midpoints(h.edges[1])[Int(uppID)]
end


##############################################################
####                RAW PLOTS + TRUTH                     ####
##############################################################
h = fit(Histogram, data, minE:deltaE:maxE)
plot(normalize(h, mode=:density), ylabel="normalized rate", xlabel="E [keV]", label="fake data")

tru_params = (muB=bkgCounts_true, sigma=sigmaTrue, muS=sigCounts_true)
plot!(x -> likelihood_a(tru_params, x) * nTot , label="true likelihood (0 signal)")


##############################################################
####                BINNED LOG LIKELIHOOD!                ####
##############################################################
binned_model = let h = h, l_a = likelihood_a, l_b = likelihood_b
    bin_counts = h.weights
    edges = h.edges[1]
    bin_centers = collect(midpoints(edges))
    bin_width = step(edges)

    n = sum(h.weights)
    logfuncdensity(function (p)
        function binned_ll(i)
            exp_counts = bin_width * l_a(p, bin_centers[i])
            logpdf(Poisson(exp_counts + eps(float(exp_counts))), bin_counts[i])
        end

        ll_a = 0.0
        for i in eachindex(bin_counts)
            ll_a += binned_ll(i)
        end

        ll_b = log(l_b(p, n))
        return ll_a + ll_b
    end)
end

##############################################################
####              UNBINNED LOG LIKELIHOOD!                ####
##############################################################
struct UnbinnedModel
    data
end

@inline DensityInterface.DensityKind(::UnbinnedModel) = IsDensity()

function DensityInterface.logdensityof(likelihood::UnbinnedModel, p)
    n = length(likelihood.data)

    ll_a = sum(log.(likelihood_a(p, likelihood.data)))
    ll_b = log(likelihood_b(p, n))
    ll_a + ll_b
end


##############################################################
####              PRIORS                                  ####
##############################################################
prior = distprod(
    muB=Uniform(0.1,45),
    muS=Uniform(1e-10, 10),
    sigma=130 # watch out sigma cannot be 0! Then Normal(Q, 0) = Inf
)

##############################################################
####              POSTERIOR  BINNED                       ####
##############################################################
binned_posterior = PosteriorMeasure(binned_model, prior)
binned_samples = bat_sample(binned_posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result

mode(binned_samples)
plot(binned_samples)

##############################################################
####              POSTERIOR  UNBINNED                     ####
##############################################################
unbinned_posterior = PosteriorMeasure(UnbinnedModel(data), prior)
unbinned_samples = bat_sample(unbinned_posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result

mode(unbinned_samples)
plot(unbinned_samples)


binned_unshaped_samples, f_flatten = bat_transform(Vector, binned_samples)
muSs = [m[2] for m in binned_unshaped_samples.v]

exp_mu_signal_90 = 
Na = 6.02214e23
m = 6.067
t = 2.5
W = 0.08192

Thalf = log(2) * (Na * m * t * 0.15 / W) / exp_mu_signal_90





posterior_struct = PosteriorMeasure(SignalBkgLikelihood(data), prior)

samples_struct = bat_sample(posterior_struct, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4))

res = samples.result

plot(res)


function credibility_interval(v::Vector{Float64}, alpha::Float64)
    # Sort the vector
    sorted_v = sort(v)
    
    # Compute the cumulative sum
    cumulative_sum = cumsum(sorted_v)
    
    # Find the index range for the desired credibility level
    index_range = searchsortedfirst(cumulative_sum, alpha)
    
    # Compute the interval based on the index range
    credibility_interval = sorted_v[1:index_range]
    
    return credibility_interval
end

credibility_interval(muSs, 0.9)

d = muSs[1:100]
sorted_mu = sort(d)
cs_mu = cumsum(sorted_mu) / sum(sorted_mu)

low = 0.0
up_idx = findfirst(x -> x >= 0.9, cs_mu)
up = sorted_mu[up_idx]

sum( sorted_mu[1:up_idx] ) / sum(sorted_mu)

for ilow = 1:length(sorted_mu), iup = 