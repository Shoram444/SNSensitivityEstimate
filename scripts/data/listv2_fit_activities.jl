using DrWatson
@quickactivate "SNSensitivityEstimate"

using Turing
using Random, Statistics, Distributions, StatsPlots
using DataFramesMeta, CSV, UnROOT, FHist
using StatsBase # For mean and std
using KernelDensity # Crucial for KDEs

# --- Data Loading and Preprocessing (same as before) ---
pint_cut = 0.02 # ns
pext_cut_low = 0.00 # ns
pext_cut_high = .5 # ns

phi_low, phi_high = 5, 180 # degrees
E_cut = 300 # keV
dy_cut = 180 # mm
dz_cut = 180 # mm

d1 = let
    f = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/phase1_50keV_Mathis_calibration/mva.root")
    d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
    @chain d begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset phi_low .< :phi .< phi_high
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :esum .> E_cut
        @subset :dy .< dy_cut
        @subset :dz .< dz_cut
        @subset :Pint .>= pint_cut
        @subset pext_cut_low .<= :Pext .<= pext_cut_high
    end
end

d1_sim = let
    f_sim = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/listv2_50keV_old_calibration/mva/bb_foil_bulk.root")
    d_sim = LazyTree(f_sim, "tree", keys(f_sim["tree"])) |> DataFrame
    @chain d_sim begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset phi_low .< :phi .< phi_high
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :esum .> E_cut
        @subset :dy .< dy_cut
        @subset :dz .< dz_cut
        @subset :Pint .>= pint_cut
        @subset pext_cut_low .<= :Pext .<= pext_cut_high
    end
end

d_K40 = let
    f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/listv2_50keV_old_calibration/mva/K40_foil_bulk.root"))
    d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
    @chain d begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset phi_low .< :phi .< phi_high
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :esum .> E_cut
        @subset :dy .< dy_cut
        @subset :dz .< dz_cut
        @subset :Pint .>= pint_cut
        @subset pext_cut_low .<= :Pext .<= pext_cut_high
    end
end

d_Pa234m = let
    f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/listv2_50keV_old_calibration/mva/Pa234m_foil_bulk.root"))
    d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
    @chain d begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset phi_low .< :phi .< phi_high
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :esum .> E_cut
        @subset :dy .< dy_cut
        @subset :dz .< dz_cut
        @subset :Pint .>= pint_cut
        @subset pext_cut_low .<= :Pext .<= pext_cut_high
    end
end

d_radon = let
    f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/listv2_50keV_old_calibration/mva/Bi214_wire_surface.root"))
    d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
    @chain d begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset phi_low .< :phi .< phi_high
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :esum .> E_cut
        @subset :dy .< dy_cut
        @subset :dz .< dz_cut
        @subset :Pint .>= pint_cut
        @subset pext_cut_low .<= :Pext .<= pext_cut_high
    end
end


# Observed data and its total count
observed_events = d1.esum
N_observed_total = length(observed_events)

# --- Define component PDFs using Kernel Density Estimation (KDE) ---
# It's good practice to set a range for KDE, especially for physical quantities
# to avoid density estimates outside meaningful bounds.
# A reasonable range would be the min/max of all your data combined, or physical limits.
min_esum = min(minimum(d1.esum), minimum(d1_sim.esum), minimum(d_K40.esum), minimum(d_Pa234m.esum), minimum(d_radon.esum))
max_esum = max(maximum(d1.esum), maximum(d1_sim.esum), maximum(d_K40.esum), maximum(d_Pa234m.esum), maximum(d_radon.esum))

# Expand the range slightly to ensure coverage and avoid boundary issues for KDEs
kde_range = range(floor(min_esum * 0.9), ceil(max_esum * 1.1), length=2048) # power of 2 for FFT in KernelDensity

# Generate KDE objects for each component
# You can optionally specify `kernel` (e.g., Normal, Epanechnikov) and `bandwidth`
# If not specified, `KernelDensity.jl` uses a default bandwidth rule (Silverman's rule).
component_kdes = [
    kde(d1_sim.esum; boundary=(min_esum, max_esum), npoints=length(kde_range)), # Signal component
    kde(d_K40.esum; boundary=(min_esum, max_esum), npoints=length(kde_range)),    # K40 component
    kde(d_Pa234m.esum; boundary=(min_esum, max_esum), npoints=length(kde_range)),# Pa234m component
    kde(d_radon.esum; boundary=(min_esum, max_esum), npoints=length(kde_range))  # Radon component
]
num_components = length(component_kdes)

# --- Define the Turing model ---
@model function unbinned_extended_model_kde(observed_events, num_components, component_kdes, N_observed_total)
    # Priors
    total_rate ~ Uniform(N_observed_total * 0.5, N_observed_total * 1.5)
    proportions ~ Dirichlet(ones(num_components))

    # Likelihood for the total number of events (Extended Likelihood part)
    N_observed_total ~ Poisson(total_rate)

    # Unbinned Likelihood for the shape of the spectrum
    for i in eachindex(observed_events)
        x = observed_events[i]
        
        # Calculate the mixture PDF value for this event x using KDEs
        mixture_pdf_val = 0.0
        for j in 1:num_components
            # pdf() is overloaded for KernelDensity.UnivariateKDE objects
            mixture_pdf_val += proportions[j] * pdf(component_kdes[j], x)
        end
        
        # Add a small epsilon to avoid log(0) issues
        Turing.@addlogprob! log(max(mixture_pdf_val, 1e-12))
    end
end

# Instantiate the model with your data
model_kde = unbinned_extended_model_kde(observed_events, num_components, component_kdes, N_observed_total)

# Sample from the posterior using NUTS
sampler = NUTS(0.65)
n_samples = 2_000
n_chains = 4
n_warmup = 1_000

chains_kde = Turing.sample(model_kde, sampler, MCMCThreads(), n_samples, n_chains; warmup_iterations = n_warmup)

# --- Analysis of results ---
println("\n--- Turing.jl Results (with KDEs) ---")
display(chains_kde)

mean_total_rate_kde = mean(chains_kde[:total_rate])
mean_proportions_kde = mean(chains_kde[:proportions])

println("\nMean total_rate (KDE): ", mean_total_rate_kde)
println("Mean proportions (KDE): ", mean_proportions_kde)

# --- Plotting the fit ---
data_hist_plot = fit(Histogram, observed_events, 0:100:3500) # Choose appropriate bins for plotting

plot(data_hist_plot; label = "Observed Data", xlabel = "Energy (keV)", ylabel = "Counts", title = "Unbinned Extended Likelihood Fit with Turing.jl (KDEs)")

# Define the fitted mixture PDF function based on the mean parameters from KDE fit
fitted_mixture_pdf_turing_kde = x -> begin
    pdf_val = 0.0
    for j in 1:num_components
        pdf_val += mean_proportions_kde[j] * pdf(component_kdes[j], x)
    end
    return pdf_val
end

plot_bin_edges = data_hist_plot.edges[1]
plot_bin_centers = [0.5 * (plot_bin_edges[i] + plot_bin_edges[i+1]) for i in 1:(length(plot_bin_edges)-1)]
plot_bin_widths = [plot_bin_edges[i+1] - plot_bin_edges[i] for i in 1:(length(plot_bin_edges)-1)]

# Calculate expected counts per bin
expected_counts_per_bin_turing_kde = [fitted_mixture_pdf_turing_kde(bc) * mean_total_rate_kde * plot_bin_widths[i]
                                      for (i, bc) in enumerate(plot_bin_centers)]

plot!(plot_bin_centers, expected_counts_per_bin_turing_kde, seriestype = :steppre, label = "Turing Fit (KDEs)", color = :red, linewidth = 2)

# Optional: Plot individual components (scaled by their fitted proportion and total rate)
# for j in 1:num_components
#     component_proportions = mean_proportions_kde[j]
#     component_fitted_pdf_kde = x -> pdf(component_kdes[j], x)
#     
#     component_expected_count_density_func = x -> component_proportions * component_fitted_pdf_kde(x) * mean_total_rate_kde
#     
#     component_expected_counts_per_bin = [component_expected_count_density_func(bc) * plot_bin_widths[i]
#                                           for (i, bc) in enumerate(plot_bin_centers)]
#     
#     plot!(plot_bin_centers, component_expected_counts_per_bin, seriestype = :steppre, label = "Component $(j) (KDE)", linewidth = 1, linestyle = :dash)
# end

# savefig(plotsdir("turing_unbinned_kde_fit.png"))