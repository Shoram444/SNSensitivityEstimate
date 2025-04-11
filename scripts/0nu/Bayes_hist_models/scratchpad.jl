using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())

using SensitivityModule
using Random, LinearAlgebra, Statistics, Distributions, Plots, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
global t0 = time()

analysisDict = Dict(
    :signal => "RH037_foil_bulk",
    :roi => (0, 175),
    :bw => 5,
    :mode => "phi"
)

ROI_a, ROI_b, bw = analysisDict[:roi][1],analysisDict[:roi][2], analysisDict[:bw]

###################################
include(scriptsdir("Params.jl"))

# Load all the processes in the directory. Function `load_processes` takes two arguments:
# 1. dir::String -> the name of the directory where the root files are stored
# 2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
all_processes = load_data_processes("fal5_8perc_Boff_TKrec_evis_bcu_J38", analysisDict[:mode], fwhm=0.0)

# declare which process is signal
signal = get_process(analysisDict[:signal], all_processes)

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
set_nTotalSim!( signal, 1e8 )
# set_nTotalSim!( signal, 1e8 )

set_signal!(background[1], false)

# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1],0.99e8)
set_nTotalSim!( background[2],1e8)
set_nTotalSim!( background[3],1e8  )
set_nTotalSim!( background[4],1e8)
set_nTotalSim!( background[5],1e8  )
set_nTotalSim!( background[6],1e8  )

@info "process initialized"
println("Processes initialized.")

for b in background
    set_bins!(b, ROI_a:bw:ROI_b+bw)
end
set_bins!(signal, ROI_a:bw:ROI_b+bw)


const bkg_hist = [(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)]  
const bkg_hist_normed = normalize.(bkg_hist, width = true)
const signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b), width = true)
const f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)

# uninformed prior for each activity
prior = NamedTupleDist(
    As = Uniform(0,1),
    Ab = [Uniform(0,1) for _ in 1:length(bkg_hist)] 
)   

t = Float64[]
while(time() - t0 < 3600*36) # do this for n hours
# for _ in 1:1 # do this for n hours
    GC.gc()
    t1 = time()
    try 
        sens = get_sens_bayes_uniform(bkg_hist, f2, signal, prior; ROI_a = ROI_a, ROI_b = ROI_b, nsteps = 10^4, nchains = 4)
        println("time to fit, t = $(time() - t1) s")
        println(sens)
        push!(t, sens)
    catch
        @warn "failed fit" 
        continue
    end
end

save_name = savename(analysisDict)

using DataFramesMeta, CSV
CSV.write("/pbs/home/m/mpetro/sps_mpetro/Projects/PhD/SNSensitivityEstimate/scripts/0nu/Bayes_hist_models/data_perBkg_10keV_binning/sensitivities_$(save_name)_$(rand(1:100000)).csv", DataFrame(thalf= t))
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

# ROI_a, ROI_b = 2000, 3400

# bkg_hist = [(restrict(b, ROI_a, ROI_b)) for b in get_bkg_counts_1D.(background)]  
# bkg_hist_normed = normalize.(bkg_hist, width = true)
# signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), ROI_a, ROI_b), width = true)

# data_hist = [get_pseudo_spectrum(b) for b in bkg_hist] |> sum 
# # data_hist = sample_histogram(bkg_hist) 

# f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)


# my_likelihood = make_hist_likelihood_uniform(data_hist, f2)

# # s_prior = 1e-5
# # # Define the concentration parameters (prior belief)
# # proportions = integral.(bkg_hist) ./ sum(integral.(bkg_hist)) .+ eps(Float64)
# # α = [s_prior, proportions...] ./ sum([s_prior, proportions...]) # Prior for each activity
# # α = ones(length(bkg_hist)+1) ./ 10#length(bkg_hist) # Prior for each activity

# # Create the Dirichlet distribution
# # # prior = NamedTupleDist(w = Dirichlet(α))
# prior = NamedTupleDist(
#     As = Uniform(0,1),
#     Ab = [
#         Uniform(0, 1),#truncated(Normal(proportions[1]), 0, 1.0),#Uniform(0, 1),
#         Uniform(0, 1),# truncated(Normal(proportions[2]), 0, 1.0)#Uniform(0, 1),
#         Uniform(0, 1),# truncated(Normal(proportions[3]), 0, 1.0)#Uniform(0, 1),
#         Uniform(0, 1),# truncated(Normal(proportions[4]), 0, 1.0)#Uniform(0, 1),
#         # truncated(Normal(proportions[5]), 0, 1.0)#Uniform(0, 1),
#         # truncated(Normal(proportions[6]), 0, 1.0)#Uniform(0, 1),
#     ]
# )

# # Apply the TransformedMeasure to wrap the Dirichlet distribution


# burnin = MCMCMultiCycleBurnin(max_ncycles = 30, nsteps_final=1000)
# mcmcalgo = MetropolisHastings(
#     weighting = RepetitionWeighting(),
#     tuning = AdaptiveMHTuning()
# )


# posterior = PosteriorMeasure(my_likelihood, prior)
# samples, _ = bat_sample(
#     posterior, 
#     MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 10^4, nchains = 4),
# )


# binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
# nDataPoints = integral(data_hist)
# muS = [par[1] for par in binned_unshaped_samples.v]

# @show exp_mu_signal_90 = quantile( muS,0.9) * nDataPoints 
# Na = 6.02214e23
# m = 6.067
# t = 2.88
# W = 0.08192
# eff= lookup(signal, ROI_a, ROI_b-50)
# Thalf = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90

# Plots.plot(samples)

# bw = 10
# xs = ROI_a:bw:ROI_b-bw/2
# # b_amps = [x for x in mean(samples)][2:end] 
# amps = mean(samples) 

# ys = [f2(amps, x) for x in xs]
# Plots.plot(data_hist, label = "data",fillrange = 1e-5, fillcolor = :match)
# Plots.plot!(xs, ys .* integral(data_hist) .* bw, linewidth = 4,  label = "fit")
# # Plots.plot!(xlims = (2000, 3200), ylims = (1e-5, 100))


# sig_amps = (As = quantile( muS,0.9), Ab = amps[2])
# ys_signal = [f2(sig_amps, x) for x in xs]
# Plots.plot!(xs, ys_signal .* integral(data_hist) .* bw, linewidth = 4, label = "signal fit")
# Plots.plot!(yscale = :log10, ylims = (1.7e-2, 1e4))
# # Plots.plot!(ylims = (0, 2))
# savefig("scripts/0nu/Bayes_hist_models/fit.svg")

# SensitivityModule.stackedhist(bkg_hist)