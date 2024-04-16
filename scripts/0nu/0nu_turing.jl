using Turing, StatsPlots, StatsBase, LinearAlgebra, DataFrames

#######################################################
#####               FUNCTIONS                   #######
#######################################################

function fit_function_exponential( Q, sigma, lambda, nSig, nBkg, x )
    (nSig * pdf( Normal(Q, sigma), x ) + nBkg * pdf( Exponential(lambda), x )) / ( nSig+nBkg )
end

function fit_function_uniform( Q, sigma, nSig, nBkg, x; minE=minE, maxE=maxE )
    (nSig * pdf( Normal(Q, sigma), x ) + nBkg * inv(maxE-minE))/(nSig+nBkg)
end

function generate_data(expected_bkg_cts_per_ROI)
    events = []
    for row in eachrow(expected_bkg_cts_per_ROI)
        bin_min = row.bins[1]
        bins_step = row.bins[2] - row.bins[1]
        num_events_in_bin = rand(Poisson(row.bExp)) # number of events in bin Poisson distributed
        e = [ bin_min + bins_step* rand() for _ in 1:num_events_in_bin ]
        append!(events, e) # each event in bin is generated uniformally through bin width
    end
    float.(events)
end

#######################################################
#####              MODEL DEFINITIONS            #######
#######################################################
@model function full_model_exponential(data; Q=Q, sigma=sigmaTrue, )
    lambda ~ Uniform(10^-2, 10.)
    mu ~ Normal(Q, sigma)
    sigma1 ~ truncated(Normal(sigmaTrue, 0.01), 0, Inf)
    # Prior distributions for the proportions
    nSignal ~ Uniform(0.0, 15.0)  # Prior for proportion of Normal data
    nBkg ~ Uniform(0.0, 100.0)  # Prior for proportion of Normal data
    
    # Likelihood of data
    for x in data
        likelihood_normal = pdf(Normal(mu, sigma1), x)
        likelihood_exponential = pdf(Exponential(lambda), x)
        likelihood_mixture = (nSignal * likelihood_normal + nBkg * likelihood_exponential) 
        likelihood_extended = pdf(Poisson(nSignal + nBkg), length(data))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end

@model function full_model_uniform(data; Q=Q, sigma=sigmaTrue, minE=minE, maxE=maxE)
    mu ~ Normal(Q, sigma)
    sigma1 ~ truncated(Normal(sigmaTrue, 0.01), 0, Inf)
    # Prior distributions for the proportions
    nSignal ~ Uniform(0.0, 15.0)  # Prior for proportion of Normal data
    nBkg ~ Uniform(0.0, 100.0)  # Prior for proportion of Normal data
    
    # Likelihood of data
    for x in data
        likelihood_normal = pdf(Normal(mu, sigma1), x)
        likelihood_uniform = inv(maxE-minE) # expectation value of uniform
        likelihood_mixture = (nSignal * likelihood_normal + nBkg * likelihood_uniform) 
        likelihood_extended = pdf(Poisson(nSignal + nBkg), length(data))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end

function generate_sample_output_uniform(expected_bkg_cts_per_ROI; α=0.9)
    data = generate_data(expected_bkg_cts_per_ROI)
    model = full_model_uniform(data)
    chains = sample(model, NUTS(1000, 0.65), MCMCThreads(), 10_000, 1)

    lambda = mean(chains[:,1, :])
    mu = mean(chains[:,2, :])
    sigma = mean(chains[:,3, :])
    nSig = mean(chains[:,4, :]) 
    nSig_90q = quantile(chains[:, 4, 1].data, 0.9)

    nBkg = mean(chains[:,5, :]) 
    nBkg_90q = quantile(chains[:, 5, 1].data, 0.9)

    ( lambda=lambda, mu=mu,sigma=sigma, nSig=nSig, 
      nSig_90q=nSig_90q, nBkg=nBkg, nBkg_90q=nBkg_90q, data=data )
end

function generate_sample_output_exponential(expected_bkg_cts_per_ROI; α=0.9)
    data = generate_data(expected_bkg_cts_per_ROI)
    model = full_model_exponential(data)
    chains = sample(model, NUTS(1000, 0.65), MCMCThreads(), 10_000, 1)

    lambda = mean(chains[:,1, :])
    mu = mean(chains[:,2, :])
    sigma = mean(chains[:,3, :])
    nSig = mean(chains[:,4, :]) 
    nSig_90q = quantile(chains[:, 4, 1], 0.9)

    nBkg = mean(chains[:,5, :]) 
    nBkg_90q = quantile(chains[:, 5, 1], 0.9)

    ( lambda=lambda, mu=mu,sigma=sigma, nSig=nSig, 
      nSig_90q=nSig_90q, nBkg=nBkg, nBkg_90q=nBkg_90q, data=data )
end

function get_tHalf(nSig_90q)
    Na = 6.02214e23
    m = 6.067
    t = 2.88
    W = 0.08192
    eff= 0.1737
    Thalf = log(2) * (Na * m * t * eff / W) / nSig_90q
end


#######################################################
#####               DATA                        #######
#######################################################
# Generate some sample data
# signal
Q = 2.99 - 2.6
sigmaTrue = 0.060
minE = 2.600 - 2.6
maxE = 3.200 - 2.6
deltaE = 0.100
expected_bkg_cts_per_ROI = DataFrame( 
    bins = [ (e1, e1+deltaE) for e1 in minE:deltaE:maxE-deltaE ],
    bExp = [  4.22099, 0.999972, 0.529423, 0.382328, 0.276759, 0.237596  ]
)

data = generate_data(expected_bkg_cts_per_ROI)

h = fit(Histogram, data, minE:deltaE:maxE)
plot(normalize(h, mode = :pdf), label = "pseudo-data", legend =:best)

#######################################################
#####              SAMPLING AND ANALYSIS SIMPLE #######
#######################################################
model = full_model_uniform(data)
chains = sample(model, NUTS(1000, 0.65), MCMCThreads(), 10_000, 4, )

mu = mean(chains[:,1,:])
sigma = mean(chains[:,2,:])
nSig = mean(chains[:,3, :]) 
nSig_90q = quantile(chains[:, 3, 1], 0.9)

nBkg = mean(chains[:,4, :]) 
nBkg_90q = quantile(chains[:, 4, 1], 0.9)

plot(chains, thickness_scaling = 1)

plot( normalize(h, mode=:pdf) )
plot!( range( minE,maxE, 100 ), x->fit_function_uniform(mu, sigma, nSig, nBkg, x  ) )


model = full_model_exponential(data )
chains = sample(model, NUTS(1000, 0.65), MCMCThreads(), 10_000, 4)

lambda = mean(chains[:,1,:])
mu = mean(chains[:,2,:])
sigma = mean(chains[:,3,:])
nSig = mean(chains[:,4, :]) 
nSig_90q = quantile(chains[:, 4, 1], 0.9)

nBkg = mean(chains[:,5, :]) 
nBkg_90q = quantile(chains[:, 5, 1], 0.9)

plot(chains, thickness_scaling = 1)

plot( normalize(h, mode=:pdf) )
plot!( range( minE, maxE, 100 ), x->fit_function_exponential(mu, sigma, lambda,nSig, nBkg, x  ) )

#######################################################
#####              SAMPLING AND ANALYSIS MULTIPLE    ##
#######################################################
lambda, mu,sigma, nSig, nSig_90q, nBkg, nBkg_90q, data1 = generate_sample_output_exponential(expected_bkg_cts_per_ROI)

T12_exponential = Float64[]
for i=1:10
    lambda, mu,sigma, nSig, nSig_90q, nBkg, nBkg_90q, data = generate_sample_output_exponential(expected_bkg_cts_per_ROI)
    push!(T12_exponential, get_tHalf(nSig_90q) )
end
mean(T12_exponential)
plot(T12_exponential)

T12_uniform = Float64[]
for i=1:10
    lambda, mu,sigma, nSig, nSig_90q, nBkg, nBkg_90q, data = generate_sample_output_uniform(expected_bkg_cts_per_ROI)
    push!(T12_uniform, get_tHalf(nSig_90q) )
end
mean(T12_uniform)
plot(T12_uniform)


@model function exp_only(data)
    lambda ~ Uniform(1e-4, 1e4)

    data .~ Exponential(lambda)
end

m = exp_only(data )
ch = sample(m, NUTS(), MCMCThreads(),10_000, 1)


l = mean(ch[:,1,:])
h1 = fit(Histogram, data , minE:deltaE:maxE)
plot(normalize(h1, mode=:pdf), label="data")
plot!(twinx(), xlims=(0,0.8),x-> pdf(Exponential(l), x), lw = 4, ls = :dash, c = 2, label = "fit", legend= :best)
