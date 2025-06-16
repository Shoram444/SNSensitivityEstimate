using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())

using SensitivityModule
using Random, LinearAlgebra, Statistics, Distributions, Plots, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
using DataFramesMeta, CSV
global t0 = time()



analysisDict = Dict(
    :signal => "bb0nuM1_foil_bulk",
    :bining => (300, 3500),
    :bin_width => 50,
    :mode => "sumE",
    :prior => 1e-3 # 1e-4 0nu, 1e-4 RH, 1e-3 M1, 1e-2 M2
)

Bin_low, Bin_high, bin_width = analysisDict[:bining][1],analysisDict[:bining][2], analysisDict[:bin_width]

###################################
include(srcdir("params/Params.jl"))

include(scriptsdir("ND/results/best_rois.jl"))

if analysisDict[:signal] == "bb0nu_foil_bulk"
    roi = bb0nu_roi
    backgrounds = [
        "bb_foil_bulk", 
        "Bi214_foil_bulk", 
        "Bi214_wire_surface", 
        "Tl208_foil_bulk", 
        "K40_foil_bulk",
        "Pa234m_foil_bulk",
        "gamma_experimental_surface"
    ]
elseif analysisDict[:signal] == "bb0nuM1_foil_bulk"
    roi = bb0nuM1_roi
        backgrounds = [
        "bb_foil_bulk", 
        "Bi214_foil_bulk", 
        "Bi214_wire_surface", 
        "Tl208_foil_bulk", 
        "K40_foil_bulk",
        "Pa234m_foil_bulk",
        "gamma_experimental_surface"
    ]
elseif analysisDict[:signal] == "bb0nuM2_foil_bulk"
    roi = bb0nuM2_roi
    backgrounds = [
        "bb_foil_bulk", 
        "Bi214_foil_bulk", 
        "Bi214_wire_surface", 
        "Tl208_foil_bulk", 
        "K40_foil_bulk",
        "Pa234m_foil_bulk",
        "gamma_experimental_surface"
    ]
else
    error("Unknown signal process: $(analysisDict[:signal])")
end

roi[:sumE] = (Bin_low, Bin_high) # update the sumE range to match the analysisDict
@show roi

# Load all the processes in the directory. Function `load_processes` takes two arguments:
# 1. dir::String -> the name of the directory where the root files are stored
# 2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
all_processes = load_data_processes("mva/fal5_TKrec_J40", analysisDict[:mode], fwhm=0.0, roi = roi)

# declare which process is signal
signal = get_process(analysisDict[:signal], all_processes) |> first

# declare background processes
background = [get_process(b, all_processes) |> first for b in backgrounds]

set_signal!(background[1], false)

@info "process initialized"
println("Processes initialized.")

for b in background
    set_bins!(b, Bin_low:bin_width:Bin_high+bin_width)
end
set_bins!(signal, Bin_low:bin_width:Bin_high+bin_width)


bkg_hist = [(restrict(b, Bin_low, Bin_high)) for b in get_bkg_counts_1D.(background)]  
bkg_hist_normed = normalize.(bkg_hist, width = true)
signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), Bin_low, Bin_high), width = true)
f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)

# uninformed prior for each activity
prior = NamedTupleDist(
    As = Uniform(1e-20, analysisDict[:prior]), # 1e-4 0nu, 1e-4 RH, 1e-3 M1, 1e-2 M2
    Ab = [Uniform(1e-20,1) for _ in 1:length(bkg_hist)] 
)   

t = Float64[]
while(time() - t0 < 3600*22) # do this for n hours
# for _ in 1:1 # do this for n hours
    GC.gc()
    t1 = time()
    try 
        sens = get_sens_bayes_uniform(bkg_hist, f2, signal, prior; ROI_a = Bin_low, ROI_b = Bin_high, nsteps = 5*10^4, nchains = 4)
        println("time to fit, t = $(time() - t1) s")
        println(sens)
        push!(t, sens)
        save_name = savename(analysisDict)
        CSV.write(scriptsdir("0nu/Bayes_hist_models/informative_priors/Cimrman_J40_ND_cuts/sensitivities_$(save_name)_$(rand(1:100000)).csv"), DataFrame(thalf= t))
    catch
        @warn "failed fit" 
        continue
    end
end


# CSV.write("scripts/0nu/Bayes_hist_models/sensitivities_$(save_name)_$(rand(1:100000)).csv", DataFrame(thalf= t))

# plot(t, st=:histogram, nbins = 10, xlabel = "Bayes sensitivity (yr)", label = "sample sensitivity")
# vline!([median(t)], label = "Median = $(round(median(t), sigdigits = 3)) yr", color = :red, linewidth = 4)
# vline!( [mean(t)], label = "Mean = $(round(mean(t), sigdigits = 3)) yr", color = :green, linewidth = 4)
# savefig("scripts/0nu/Bayes_hist_models/histo_sens.svg")

# plot(samples)
# # plot(signal_hist_normed)

# median(t)

# scatter(t)


# for b in background
#     set_bins!(b, 0:10:3500)
# end
# set_bins!(signal, 0:10:3500)

# Bin_low, Bin_high = 2000, 3400

# bkg_hist = [(restrict(b, Bin_low, Bin_high)) for b in get_bkg_counts_1D.(background)]  
# bkg_hist_normed = normalize.(bkg_hist, width = true)
# signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), Bin_low, Bin_high), width = true)

# data_hist = [get_pseudo_spectrum(b) for b in bkg_hist] |> sum 
# data_hist = sample_histogram(bkg_hist) 

# f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)


# my_likelihood = make_hist_likelihood_uniform(data_hist, f2)

# # # s_prior = 1e-5
# # # # Define the concentration parameters (prior belief)
# # # proportions = integral.(bkg_hist) ./ sum(integral.(bkg_hist)) .+ eps(Float64)
# # # α = [s_prior, proportions...] ./ sum([s_prior, proportions...]) # Prior for each activity
# # # α = ones(length(bkg_hist)+1) ./ 10#length(bkg_hist) # Prior for each activity

# # # Create the Dirichlet distribution
# # # # prior = NamedTupleDist(w = Dirichlet(α))
# # prior = NamedTupleDist(
# #     As = Uniform(0,1),
# #     Ab = [
# #         Uniform(0, 1),#truncated(Normal(proportions[1]), 0, 1.0),#Uniform(0, 1),
# #         Uniform(0, 1),# truncated(Normal(proportions[2]), 0, 1.0)#Uniform(0, 1),
# #         Uniform(0, 1),# truncated(Normal(proportions[3]), 0, 1.0)#Uniform(0, 1),
# #         Uniform(0, 1),# truncated(Normal(proportions[4]), 0, 1.0)#Uniform(0, 1),
# #         # truncated(Normal(proportions[5]), 0, 1.0)#Uniform(0, 1),
# #         # truncated(Normal(proportions[6]), 0, 1.0)#Uniform(0, 1),
# #     ]
# # )

# # # Apply the TransformedMeasure to wrap the Dirichlet distribution


# burnin = MCMCMultiCycleBurnin(max_ncycles = 30, nsteps_final=1000)
# mcmcalgo = MetropolisHastings(
#     weighting = RepetitionWeighting(),
#     tuning = AdaptiveMHTuning()
# )


# posterior = PosteriorMeasure(my_likelihood, prior)
# samples, _ = bat_sample(
#     posterior, 
#     MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 5*10^4, nchains = 4),
# )


# binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
# nDataPoints = integral(data_hist)
# muS = [par[1] for par in binned_unshaped_samples.v]

# @show exp_mu_signal_90 = quantile( muS,0.9) * nDataPoints 
# Na = 6.02214e23
# m = 6.067
# t = 2.88
# W = 0.08192
# eff= lookup(signal, Bin_low, Bin_high-50)
# Thalf = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90

# Plots.plot(samples, size = (2000, 1600))
# Plots.plot(samples, :As)
# savefig("notebooks/Sensitivity_exotic_Bayes_nDim/samples_As.png")

# bin_width = 1#analysisDict[:bin_width]
# xs = Bin_low:bin_width:Bin_high-bin_width/2
# # b_amps = [x for x in mean(samples)][2:end] 
# amps = mean(samples) 

# ys = [f2(amps, x) for x in xs]
# Plots.plot(data_hist, label = "data",fillrange = 1e-5, fillcolor = :match)
# Plots.plot!(xs, ys .* integral(data_hist) .* bin_width * 5, linewidth = 4,  label = "fit")
# # # Plots.plot!(xlims = (2000, 3200), ylims = (1e-5, 100))


# sig_amps = (As = quantile( muS,0.9), Ab = amps[2])
# ys_signal = [f2(sig_amps, x) for x in xs]
# Plots.plot!(xs, ys_signal .* integral(data_hist) .* bin_width *5, linewidth = 4, label = "signal fit")
# Plots.plot!(yscale = :log10, ylims = (1.7e-2, 1e4))
# Plots.plot!(ylims = (0, 2))
# # savefig("scripts/0nu/Bayes_hist_models/fit.svg")

# # SensitivityModule.stackedhist(bkg_hist)