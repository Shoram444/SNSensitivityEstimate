using DrWatson
@quickactivate "SNSensitivityEstimate"

using BAT
using Random, Statistics, Distributions, StatsPlots
using DataFramesMeta, CSV, UnROOT, FHist
using StatsBase # For mean and std
using BinnedModels, ValueShapes
using MeasureBase # Crucial for KDEs

bin_low, bin_width, bin_high = 0, 50, 3500

# --- Data Loading and Preprocessing (same as before) ---
pint_cut = 0.02 # ns
pext_cut_low = 0.00 # ns
pext_cut_high = .5 # ns

phi_low, phi_high = 5, 180 # degrees
E_cut = 300 # keV
dy_cut = 180 # mm
dz_cut = 180 # mm

d1 = let
    f = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v2_phase1_50keV_Xalbat_calib/mva.root")
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
    f_sim = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/bb_foil_bulk.root")
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
    f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/K40_foil_bulk.root"))
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
    f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/Pa234m_foil_bulk.root"))
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
    f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/Bi214_wire_surface.root"))
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


data_hist = fit(Histogram, d1.esum, bin_low:bin_width:bin_high)

hist_sim = let
    hists = FHist.Hist1D[]
    push!(hists, normalize(Hist1D(d1_sim.esum; binedges= bin_low:bin_width:bin_high)))
    push!(hists, normalize(Hist1D(d_K40.esum; binedges= bin_low:bin_width:bin_high)))
    push!(hists, normalize(Hist1D(d_Pa234m.esum; binedges= bin_low:bin_width:bin_high)))
    push!(hists, normalize(Hist1D(d_radon.esum; binedges= bin_low:bin_width:bin_high)))
    hists
end


# Define the total observed events for the Poisson term
N_observed_total = sum(data_hist.weights)


function f_expectation(par, hist_sim, epsilon_count = 1e-9) # A very small positive number
    proportions = par.rate
    num_sources = length(proportions)

    @assert num_sources == length(hist_sim) "Number of proportions must match number of templates"

    return x -> begin
        pdf_value = 0.0
        for i in 1:num_sources
            pdf_value += proportions[i] * lookup(hist_sim[i], x)
        end
        # Add epsilon to prevent zero PDF values
        return max(pdf_value, epsilon_count)
    end
end

model = binned_model(par -> f_expectation(par, hist_sim), data_hist.edges)
likelihood = Likelihood(model, data_hist.weights)

# prior = NamedTupleDist((
#     rate = Dirichlet(ones(length(hist_sim))),  # Dirichlet over the N-simplex
# ))

prior = distprod(
    rate = Dirichlet(ones(length(hist_sim))),  # Dirichlet over the N-simplex
)


posterior = PosteriorMeasure(likelihood, prior)

results = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 5*10^4, nchains = 6))
samples = results.result

println("Mode: $(mode(samples))")
println("Mean: $(mean(samples))")
println("Median: $(median(samples))")
println("Stddev: $(std(samples))")

unshaped_samples, f_flatten = bat_transform(Vector, samples)

par_cov = cov(unshaped_samples)
println("Covariance: $par_cov")

plot(samples)

function fit_function(p::NamedTuple{(:rate,), T}, x::Real, hist_sim::Vector{Hist1D}) where T
    sum(p.rate[i] * lookup(hist_sim[i], x) for i in eachindex(hist_sim))
end

fit_function(p::NamedTuple{(:rate,), T}, x::Real) where T = fit_function(p, x, hist_sim)

mean_vals = mean(samples)

plot(normalize(data_hist) ; label = "Data", xlabel = "Energy (keV)", ylabel = "Counts", title = "Fit to Data")
plot!(0:10:3490, x -> fit_function(mean_vals, x), label = "Fit", color = :red, linewidth = 4)
# plot(0:10:3490, fit_function, samples)


# Reshape and collect samples
shaped_samples = [x.v for x in collect(samples)]

plot(0:10:3490, x -> fit_function(shaped_samples[1], x), label="Sample 1")

for i in rand(1:length(shaped_samples), 50)
    plot!(0:10:3490, x -> fit_function(shaped_samples[i], x), label="", alpha=0.2, color=:gray)
end
current()





# Evaluate fit_function5 for many samples and many x
xs = 0:10:3490
n_x = length(xs)
n_samples = length(shaped_samples)

# Collect all predictions for each x
ys = [fit_function(p, x) for p in shaped_samples, x in xs]  # size: (n_samples, n_x)

# Transpose to get ys_per_x[i] = all sample predictions at x = xs[i]
ys_per_x = eachcol(ys)

# Compute stats for each x
med = [median(y) for y in ys_per_x]
low_1σ = [quantile(y, 0.16) for y in ys_per_x]
high_1σ = [quantile(y, 0.84) for y in ys_per_x]
low_2σ = [quantile(y, 0.025) for y in ys_per_x]
high_2σ = [quantile(y, 0.975) for y in ys_per_x]

# Plot!
plot(xs, med; label="Median fit", lw=2, color=:red)
plot!(xs, low_1σ, ribbon=(high_1σ .- low_1σ), fillalpha=0.4, label="68% band", color=:red)
plot!(xs, low_2σ, ribbon=(high_2σ .- low_2σ), fillalpha=0.2, label="95% band", color=:red)

# Optional: overlay normalized data histogram
norm_data_hist = normalize(data_hist)
bar!(norm_data_hist.edges[1][1:end-1], norm_data_hist.weights;
     label="Data", legend=:topright, alpha=0.5, bar_width=100)

current()