using Turing, StatsPlots, StatsBase, LinearAlgebra

# Generate some sample data
# signal
mu1_true = 5
sigma1_true = 0.5
n1_true = 0

# background
lambda_true = 5
n2_true = 25


data = vcat(rand(Normal(mu1_true, sigma1_true), n1_true),
             rand(Exponential(lambda_true), n2_true))

h = fit(Histogram, data, nbins = 10)
plot(normalize(h, mode = :density), label = "data")

# Run the sampler
@model function full_model(data)
    mu1 ~ truncated(Normal(5, 0.5), 1e-6, 20)
    sigma1 ~ truncated(Normal(0.5, 0.1), 0.1, 2.0)
    lambda ~ truncated(Normal(5, 2), 1.0, 100.0)

    # Prior distributions for the proportions
    p_normal ~ Uniform(0, length(data))  # Prior for proportion of Normal data
    p_exponential ~ Uniform(0, length(data))  # Prior for proportion of Normal data
    # p_exponential = 1 - p_normal  # Calculate proportion of Exponential data
    
    # Likelihood of data
    for x in data
        likelihood_normal = pdf(Normal(mu1, sigma1), x)
        likelihood_exponential = pdf(Exponential(lambda), x)
        likelihood_mixture = (p_normal * likelihood_normal + p_exponential * likelihood_exponential) / (p_normal + p_exponential)
        likelihood_extended = pdf(Poisson(p_normal + p_exponential), length(data))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end

model = full_model(data)
chains3 = sample(model, NUTS(), MCMCThreads(), 10_000, 4)

function fit_function( mu1, sigma1, lambda, n1, n2, x )
    n1 * pdf( Normal(mu1, sigma1), x ) + n2 * pdf( Exponential(lambda), x )
end

mu1 = mean(chains3[:,1, :])
sigma1 = mean(chains3[:,2, :])
lambda = mean(chains3[:,3, :])
n1 = mean(chains3[:,4, :]) 
n1_90quantile = quantile(chains3[:, 4, :], [0., 0.9])

n2 = mean(chains3[:,5, :]) 
n2_90quantile = quantile(chains3[:, 5, :], [0., 0.9])


h = fit(Histogram, data, nbins = 10)

plot( h , label ="data" , legend= :best)
plot!( x -> fit_function( mu1, sigma1, lambda, n1, n2, x ), label ="fit" )


current()

plot(chains3, thickness_scaling = 0.8)
