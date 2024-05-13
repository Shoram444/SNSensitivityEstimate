using DrWatson
@quickactivate 

using BAT, Distributions, StatsPlots, DataFrames, StatsBase, LinearAlgebra
using DensityInterface, IntervalSets, BinnedModels

Q = 2.99
sigmaTrue = 0.060
minE = 2.700
maxE = 3.500
deltaE = 0.100
expected_bkg_cts_per_ROI = DataFrame( 
    bins = [ (e1, e1+deltaE) for e1 in minE:deltaE:maxE-deltaE ],
    bExp = [  0.999972, 0.529423, 0.382328, 0.276759, 0.237596, 0.207903, 0.172914, 0.125756  ] 
)

data = let b=expected_bkg_cts_per_ROI
    events = []
    for row in eachrow(b)
        bin_min = row.bins[1]
        bins_step = row.bins[2] - row.bins[1]
        num_events_in_bin = rand(Poisson(row.bExp)) # number of events in bin Poisson distributed
        @show e = [ bin_min + bins_step* rand() for _ in 1:num_events_in_bin ]
        append!(events, e) # each event in bin is generated uniformally through bin width
    end
    events
end

data = float.(data)
h = fit(Histogram, data, minE:deltaE:maxE)
plot(h)

# function likelihood_a(p::NamedTuple{(:muB, :muS, :sigma)}, x::Real; minE=minE, maxE=maxE, Q=Q)
#     (p.muB * pdf(Uniform(minE, maxE), x)  +
#      p.muS * pdf(Normal(Q, p.sigma),  x)) *
#     inv( p.muB + p.muS ) 
# end


function likelihood_a(muB::Real, muS::Real, sigma::Real, lambda::Real, x::Real; minE=minE, maxE=maxE, Q=Q)
    (muB * pdf(Exponential(lambda), x)  +
     muS * pdf(Normal(Q, sigma),  x)) *
    inv( muB + muS ) 
end

function likelihood_a(p::NamedTuple{(:muB, :muS, :sigma, :lambda)}, x::Real; minE=minE, maxE=maxE, Q=Q)
    likelihood_a(p.muB, p.muS, p.sigma, p.lambda, x; minE=minE, maxE=maxE, Q=Q)
end

function likelihood_b(p, n)
    pdf(Poisson(p.muB + p.muS), n)
end

function fit_function(p, x; minE=minE, maxE=maxE, Q=Q, deltaE=deltaE)
    deltaE*(p.muB * pdf( Exponential(p.lambda), x ) + p.muS * pdf( Normal(Q, p.sigma), x) )
end

"""
    Stupid integration for getting upper limit - requires monotonely decreasing samples in vector form!
"""
function get_interval_upper(data, CL; nbins=100)
    h = fit(Histogram, data; nbins = nbins)
    cs = cumsum(h.weights) ./ sum(h.weights)

    uppID = findfirst(x -> x >= CL, cs)
    midpoints(h.edges[1])[Int(uppID)]
end


##############################################################
####                RAW PLOTS + TRUTH                     ####
##############################################################
plot(h, ylabel="counts", xlabel="E [MeV]", label="fake data", legend =:best)
plot(normalize(h, mode=:pdf), ylabel="normalized rate", xlabel="E [MeV]", label="fake data", legend =:best)



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
    muB, muS, sigma, lambda = p

    ll_a = sum(log.(likelihood_a.(muB, muS, sigma, lambda, likelihood.data)))
    ll_b = log(likelihood_b(p, n))
    ll_a + ll_b
end


##############################################################
####              PRIORS                                  ####
##############################################################
priorBkg = distprod(
    muB=Uniform(1e-5, 1e3),
    muS=0.0,
    sigma=sigmaTrue, # watch out sigma cannot be 0! Then Normal(Q, 0) = Inf
    lambda = Uniform(1e-5,1e2 )
)

prior = distprod(
    muB=Uniform(1e-5, 1e4), #Uniform(1e-5, 45),
    muS=Uniform(1e-5, 1e4), #Uniform(1e-5, 15),
    sigma=sigmaTrue, # watch out sigma cannot be 0! Then Normal(Q, 0) = Inf
    lambda = Uniform(1e-5,1e5 )
)

##############################################################
####              POSTERIOR  BINNED BKG ONLY              ####
##############################################################
binned_posterior_Bkg = PosteriorMeasure(binned_model, priorBkg)
binned_samples_Bkg = bat_sample(binned_posterior_Bkg, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result

mode(binned_samples_Bkg)
plot(binned_samples_Bkg)

##############################################################
####              POSTERIOR  BINNED                       ####
##############################################################
binned_posterior = PosteriorMeasure(binned_model, prior)
binned_samples = bat_sample(binned_posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=6)).result

mode(binned_samples)
plot(binned_samples)

binned_unshaped_samples, f_flatten = bat_transform(Vector, binned_samples)
muSsBinned = [m[2] for m in binned_unshaped_samples.v]


exp_mu_signal_90 = get_interval_upper(muSsBinned, 0.9)
Na = 6.02214e23
m = 6.067
t = 2.88
W = 0.08192
eff= 0.1737
Thalf = log(2) * (Na * m * t * 0.15 / W) / exp_mu_signal_90

binned_params=mode(binned_samples)

plot(binned_samples, size = (2000, 1600), thickness_scaling=1.6)
safesave(plotsdir("Bayessian0nu", "binned_samples_BAT.png" ), current())
safesave(plotsdir("Bayessian0nu", "binned_samples_BAT.pdf" ), current())

let f= fit_function
    plot(h, bins=(minE:deltaE:maxE),legend= :best, label= "data", widen =:false, xlims =(minE, maxE))
    plot!(
        range(minE, maxE, length=500), 
        f, 
        binned_samples, 
        xlabel = "Energy [keV]", 
        ylabel = "Background distribution", safesave(plotsdir("Bayessian0nu", "binned_fit.png" ), current())
        # sa
        title= "binned_fit: Exponential + Normal likelihood ",
        colors = [2,3,4], 
        fa = 0.6, 
        interval_labels = ["1σ", "2σ", "3σ"],
        size = (1200, 800),
        thickness_scaling= 1.6,
        legend=:best,
        dpi = 200
    )
    # safesave(plotsdir("Bayessian0nu", "binned_fit.png" ), current())
    # safesave(plotsdir("Bayessian0nu", "binned_fit.pdf" ), current())
    current()
end

##############################################################
####              POSTERIOR  UNBINNED                     ####
##############################################################
unbinned_posterior = PosteriorMeasure(UnbinnedModel(data), prior)
unbinned_samples = bat_sample(unbinned_posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4)).result

mode(unbinned_samples)
plot(unbinned_samples, size = (2000, 1600), thickness_scaling=1.6)
safesave(plotsdir("Bayessian0nu", "unbinned_samples.png" ), current())
safesave(plotsdir("Bayessian0nu", "unbinned_samples.pdf" ), current())

unbinned_unshaped_samples, f_flatten = bat_transform(Vector, unbinned_samples)
muSsUnbinned = [m[2] for m in unbinned_unshaped_samples.v]


exp_mu_signal_90 = get_interval_upper(muSsUnbinned, 0.9)
Na = 6.02214e23
m = 6.067
t = 2.88
W = 0.08192
eff= 0.1737
Thalf = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90

unbinned_params=mode(unbinned_samples)

let f= fit_function(p, x) = likelihood_a(p,x)
    plot(normalize(h, mode =:pdf), legend= :best, label= "data", widen =:false, xlims =(minE, maxE))
    plot!(
        range(minE, maxE, length=500), 
        f, 
        unbinned_samples, 
        xlabel = "Energy [keV]", 
        ylabel = "Background distribution", 
        title= "unbinned_fit: Exponential + Normal likelihood ",
        colors = [2,3,4], 
        fa = 0.6, 
        interval_labels = ["1σ", "2σ", "3σ"],
        size = (1200, 800),
        thickness_scaling= 1.6,
        legend=:best,
        dpi = 200
    )
    safesave(plotsdir("Bayessian0nu", "unbinned_fit.png" ), current())
    safesave(plotsdir("Bayessian0nu", "unbinned_fit.pdf" ), current())
    current()
end

ll = inv(mode(unbinned_samples).lambda)
plot!( 2.5:0.1:3.5, x-> ll*exp(-ll*x) )

normalize(h, mode=:pdf) 