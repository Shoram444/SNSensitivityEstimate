using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())

using SensitivityModule
using Random, LinearAlgebra, Statistics, Distributions, Plots, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
global t0 = time()


my_pdf(h::Hist1D, x::Real) = ismissing(lookup(h, x)) ? 0.0 : lookup(h, x)


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
            expected_counts = bin_width * f(Vector(p.w), bin_center) * n
            expected_counts == 0.0 && continue

            log_pdf_poisson(expected_counts, observed_counts) == -Inf && continue
            # Add log of Poisson probability for bin:
            ll_value += log_pdf_poisson(expected_counts, observed_counts)
            
        end
        ll_value
        return ll_value
    end )
end

###################################
include(scriptsdir("Params.jl"))

# Load all the processes in the directory. Function `load_processes` takes two arguments:
# 1. dir::String -> the name of the directory where the root files are stored
# 2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
all_processes = load_data_processes("fal5_8perc_Boff_TKrec_evis_bcu_J38", "sumE", fwhm=0.0)

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
set_nTotalSim!( background[1],1e8)
set_nTotalSim!( background[2],1e8)
set_nTotalSim!( background[3],1e8  )
set_nTotalSim!( background[4],1e8)
set_nTotalSim!( background[5],1e8  )
set_nTotalSim!( background[6],1e8  )

@info "process initialized"
println("Processes initialized.")

for b in background
    set_bins!(b, 0:10:3500)
end
set_bins!(signal, 0:10:3500)

bkg_hist = [b for b in get_bkg_counts_1D.(background)] 
bkg_hist_normed = normalize.(bkg_hist, width = true)
signal_hist_normed = normalize(get_bkg_counts_1D(signal), width = true)

function f(pars::Vector{Float64}, x::Real, s_hist::Hist1D, b_hists::Vector{<:Hist1D})
    sig = pars[1]*s_hist
    hists = [ pars[i]*b_hists[i-1] for i in 2:length(pars)]
    th = normalize(
        sum(vcat(sig, hists)), 
        width = true
    )
    
    return my_pdf(th, x) 
end
f1(pars::Vector{Float64}, x::Real) = f(pars, x, signal_hist_normed, bkg_hist_normed)


for b in background
    set_bins!(b, 0:50:3500)
end
set_bins!(signal, 0:50:3500)

function sample_histogram(h::Hist1D)
    data_bkg = [first(FHist.sample(h)) for i=1:rand(Poisson(round(Int, integral(h))))] 
    Hist1D( data_bkg; binedges= binedges(h) )
end

function get_sens_bayes(background::Vector{<:DataProcess}, signal::DataProcess)
    ROI_a, ROI_b = 1000, 3400

    bkg_hist = [(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)] 
    bkg_hist_normed = normalize.(bkg_hist, width = true)
    signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b), width = true)
    
    data_hist = [sample_histogram(b) for b in bkg_hist] |> sum 

    my_likelihood = make_hist_likelihood(data_hist, f1)

    # uninformed prior for each activity
    α = ones(length(bkg_hist)) ./ length(bkg_hist) 
    # Create the Dirichlet distribution
    prior = NamedTupleDist(w = Dirichlet(α))
    
    burnin = MCMCMultiCycleBurnin(nsteps_final=1000)
    mcmcalgo = MetropolisHastings(
        weighting = RepetitionWeighting(),
        tuning = AdaptiveMHTuning()
    )

    
    posterior = PosteriorMeasure(my_likelihood, prior)
    samples, _ = bat_sample(posterior, MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 10^4, nchains = 4, strict = false))


    binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
    nDataPoints = integral(data_hist)
    muS = [par[1] * nDataPoints for par in binned_unshaped_samples.v]

    @show exp_mu_signal_90 = quantile( muS,0.9) 
    Na = 6.02214e23
    m = 6.067
    t = 2.88
    W = 0.08192
    eff= lookup(signal, ROI_a, ROI_b)
    return log(2) * (Na * m * t * eff / W) / exp_mu_signal_90
end

t = Float64[]
# while(time() - t0 < 3600*10) # do this for n hours
for _ in 1:1 # do this for n hours
    GC.gc()
    println("elapsed time = $(time() - t0) s")

    sens = get_sens_bayes(background, signal)
    println(sens)
    push!(t, sens)
end

using DataFramesMeta, CSV
CSV.write("/pbs/home/m/mpetro/sps_mpetro/Projects/PhD/SNSensitivityEstimate/scripts/0nu/Bayes_hist_models/data_perBkg_10keV_binning/sensitivities_$(rand(1:100000)).csv", DataFrame(thalf= t))
# CSV.write("scripts/0nu/Bayes_hist_models/sensitivities_$(rand(1:100000)).csv", DataFrame(thalf= t))

# plot(t, st=:histogram, nbins = 10, xlabel = "Bayes sensitivity (yr)", label = "sample sensitivity")
# vline!([median(t)], label = "Median = $(round(median(t), sigdigits = 3)) yr", color = :red, linewidth = 4)
# vline!( [mean(t)], label = "Mean = $(round(mean(t), sigdigits = 3)) yr", color = :green, linewidth = 4)
# savefig("scripts/0nu/Bayes_hist_models/histo_sens.svg")

# plot(samples)
# # plot(signal_hist_normed)

# median(t)

# scatter(t)



# ROI_a, ROI_b = 1000, 3400

# bkg_hist = [(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)] 
# bkg_hist_normed = normalize.(bkg_hist, width = true)
# signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b), width = true)

# data_hist = [sample_histogram(b) for b in bkg_hist] |> sum 

# my_likelihood = make_hist_likelihood(data_hist, f1)

# # s_prior = 1e-5
# # Define the concentration parameters (prior belief)
# # proportions = integral.(bkg_hist) ./ sum(integral.(bkg_hist)) .+ eps(Float64)
# # α = [s_prior, proportions...] ./ sum([s_prior, proportions...]) # Prior for each activity
# α = ones(length(bkg_hist)) ./ length(bkg_hist) # Prior for each activity

# # Create the Dirichlet distribution
# prior = NamedTupleDist(w = Dirichlet(α))

# # Apply the TransformedMeasure to wrap the Dirichlet distribution


# burnin = MCMCMultiCycleBurnin(nsteps_final=1000)
# mcmcalgo = MetropolisHastings(
#     weighting = RepetitionWeighting(),
#     tuning = AdaptiveMHTuning()
# )


# posterior = PosteriorMeasure(my_likelihood, prior)
# samples, _ = bat_sample(posterior, MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 10^4, nchains = 4, strict = false))


# binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
# nDataPoints = integral(data_hist)
# muS = [par[1] * nDataPoints for par in binned_unshaped_samples.v]

# @show exp_mu_signal_90 = quantile( muS,0.9) 
# Na = 6.02214e23
# m = 6.067
# t = 2.88
# W = 0.08192
# eff= lookup(signal, ROI_a, ROI_b)
# Thalf = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90

# plot(samples)

# # bw = 50
# # xs = 0:bw:3350
# # b_amps = [x for x in mean(samples).w][2:end]
# # ys = [f1((vcat(exp_mu_signal_90/nDataPoints, b_amps)), x) for x in xs]
# # plot(data_hist, label = "data",fillrange = 1e-5, fillcolor = :match)
# # plot!(xs, ys .* integral(data_hist) .* bw, linewidth = 4, ylims = (0, 1e4), label = "fit")
# # plot!(xlims = (2000, 3200), ylims = (1e-5, 100))


# # ys_signal = bincounts(signal_hist_normed) .* bw .* exp_mu_signal_90
# # plot!(xs, ys_signal, linewidth = 4, label = "signal")
# # plot!(yscale = :log10, ylims = (1.7e-1, 1e4))
# # savefig("scripts/0nu/Bayes_hist_models/fit.svg")