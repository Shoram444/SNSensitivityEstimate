using Turing
using StatsPlots
using Distributions
using Random
using ReverseDiff

# Set seed for reproducibility
Random.seed!(42)

# Generate synthetic data for the example
# True distribution parameters
true_mu = 50
true_sigma = 10
n_samples = 500

# Generate true energies
true_energies = rand(Normal(true_mu, true_sigma), n_samples)

# Detector response: add Gaussian noise
detector_sigma = 15
observed_energies = true_energies .+ rand(Normal(0, detector_sigma), n_samples)

# Plot the true and observed energy distributions
histogram(true_energies, bins=0:5:100, alpha=0.5, label="True Energies", xlabel="Energy", ylabel="Frequency")
histogram!(observed_energies, bins=0:5:100, alpha=0.5, label="Observed Energies")
# Bayesian Unfolding using Turing
@model function unfolding_model(observed_energies, detector_sigma)
    # Priors for the true distribution parameters
    mu ~ Normal(0, 100)
    sigma ~ truncated(Normal(0, 10), 0, Inf)

    # Prior for true energies
    true_energy_dist = Vector{Real}(undef, length(observed_energies))
    for i in 1:length(observed_energies)
        true_energy_dist[i] ~ Normal(mu, sigma)
    end

    # Likelihood of the observed data
    for i in 1:length(observed_energies)
        observed_energies[i] ~ Normal(true_energy_dist[i], detector_sigma)
    end
end

# Instantiate the model
model = unfolding_model(observed_energies, detector_sigma)
Turing.setadbackend(:reversediff)
# Sample from the posterior distribution
chain = sample(model, NUTS(0.9; max_depth = 12), 2000)

# Extract the posterior samples for the true energies
true_energy_samples = Array(chain)[:, 3:end]

# Plot the posterior distribution
plot(1, c=:grey,label="unfolded energies", lw= 6)
for i = 1800:2000
    density!(true_energy_samples[i,:], alpha=0.2, label="", c=:grey)
end
density!(true_energies,  alpha=1, label="true Energies", c= 1, lw=6)
density!(observed_energies,  alpha=1, label="observed Energies", c= :black, lw=6)

