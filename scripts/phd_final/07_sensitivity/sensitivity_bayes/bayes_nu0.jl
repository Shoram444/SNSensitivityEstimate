using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
using DataFramesMeta, CSV
using Plots
global t0 = time()



analysisDict = Dict(
    :signal => "bb0nuM1_foil_bulk", #"%SIGNAL",
    :bining => (300, 3500),#(%BINLOW, %BINHIGH),
    :bin_width => 100, #%BINWIDTH,
    :mode => "sumE", #"%MODE",
    :prior => 1e-3,#%PRIOR # 1e-4 0nu, 1e-4 RH, 1e-3 M1, 1e-2 M2
    :radon_tag => 1
)

Bin_low, Bin_high, bin_width = analysisDict[:bining][1],analysisDict[:bining][2], analysisDict[:bin_width]

###################################
include(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/results/20260123_results/best_rois.jl"))

if analysisDict[:signal] == "bb0nu_foil_bulk" && analysisDict[:radon_tag] == 1
    roi = bb0nu_roi_radon1

elseif analysisDict[:signal] == "bb0nu_foil_bulk" && analysisDict[:radon_tag] == 2
    roi = bb0nu_roi_radon2
    
elseif analysisDict[:signal] == "bb0nuM1_foil_bulk" && analysisDict[:radon_tag] == 1
    roi = bb0nuM1_roi_radon1

elseif analysisDict[:signal] == "bb0nuM2_foil_bulk" && analysisDict[:radon_tag] == 1
    roi = bb0nuM2_roi_radon1
else
    error("Unknown signal process: $(analysisDict[:signal])")
end

backgrounds = [    
    # 2nubb
    "bb_foil_bulk", 

    # internal contamination
    "Bi214_foil_bulk", 
    "K40_foil_bulk", 
    "Pa234m_foil_bulk", 
    "Tl208_foil_bulk", 

    # radon
    "Bi210_wire_surface",     
    "Bi214_wire_surface", 

    # detector contamination
    "Bi214_PMT_glass_bulk",
    "Bi214_PMT5_glass_bulk",
    "K40_PMT_glass_bulk", 
    "K40_PMT5_glass_bulk", 
    "Tl208_PMT_glass_bulk", 
    "Tl208_PMT5_glass_bulk", 
    "Bi214_calo_8inch_scin_bulk", 
    "K40_calo_8inch_scin_bulk", 

    # gamma
    "Tl208_hall_bulk", 
    "Bi214_hall_bulk", 
    "K40_hall_bulk",
]

roi[:sumE] = (Bin_low, Bin_high) # update the sumE range to match the analysisDict
@show roi

# Load all the processes in the directory. Function `load_processes` takes two arguments:
# 1. dir::String -> the name of the directory where the root files are stored
# 2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
all_processes = load_data_processes(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), analysisDict[:mode], fwhm=0.0, roi = roi)

neutron_data_dir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/neutrons_jan_2026/")

include(joinpath(neutron_data_dir, "read_neutrons_1D.jl"))
neutron_processes = load_neutron_process1D(neutron_data_dir, analysisDict[:mode]; roi=roi)
for p in neutron_processes
    set_activity!(p, p.activity / 3)
    println("Neutron process: ", p.isotopeName, " nTotalSim: ", p.nTotalSim, " activity: ", p.activity)
end
background = vcat(backgrounds, neutron_processes)

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


function get_combined_bkg_hists(processes, Bin_low, Bin_high, bin_width)
    get_b_restricted(b) = restrict(get_bkg_counts_1D(b), Bin_low, Bin_high+bin_width)
    bins = Bin_low:bin_width:Bin_high+bin_width

    h_bb = Hist1D(;binedges = bins)
    h_internal = Hist1D(;binedges = bins)
    h_radon = Hist1D(;binedges = bins)
    h_Bi210 = Hist1D(;binedges = bins)
    h_detector = Hist1D(;binedges = bins)
    h_external = Hist1D(;binedges = bins)

    bb = ["bb_foil_bulk"]
    internals = ["Bi214_foil_bulk", "K40_foil_bulk", "Pa234m_foil_bulk", "Tl208_foil_bulk"]
    radon = ["Bi214_wire_surface"]
    Bi210 = ["Bi210_wire_surface"]
    detectors = ["Bi214_PMT_glass_bulk","Bi214_PMT5_glass_bulk","K40_PMT_glass_bulk","K40_PMT5_glass_bulk","Tl208_PMT_glass_bulk","Tl208_PMT5_glass_bulk","Bi214_calo_8inch_scin_bulk","K40_calo_8inch_scin_bulk", ]

    for p in processes
        h = get_b_restricted(p)

        if p.isotopeName in internals && integral(h) >0
            h_internal = h_internal + h
        elseif p.isotopeName in bb && integral(h) >0
            h_bb = h_bb + h
        elseif p.isotopeName in radon && integral(h) >0
            h_radon = h_radon + h
        elseif p.isotopeName in Bi210 && integral(h) >0
            h_Bi210 = h_Bi210 + h
        elseif p.isotopeName in detectors && integral(h) >0
            h_detector = h_detector + h
        elseif integral(h) >0
            h_external = h_external + h
        end
    end
    return vcat(h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external)
end

bkg_hist = get_combined_bkg_hists(background, Bin_low, Bin_high, bin_width)


bkg_hist_normed = normalize.(bkg_hist, width = true)
signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), Bin_low, Bin_high+bin_width), width = true)
f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)

# uninformed prior for each activity
prior = NamedTupleDist(
    As = Uniform(1e-20, analysisDict[:prior]), # 1e-4 0nu, 1e-4 RH, 1e-3 M1, 1e-2 M2
    Ab = [Uniform(1e-20,1) for _ in 1:length(bkg_hist)] 
)   

t = Float64[]
# while(time() - t0 < 3600*22) # do this for n hours
# # for _ in 1:1 # do this for n hours
#     GC.gc()
#     t1 = time()
#     try 
#         sens = get_sens_bayes_uniform(bkg_hist, f2, signal, prior; ROI_a = Bin_low, ROI_b = Bin_high, nsteps = 5*10^4, nchains = 4)
#         println("time to fit, t = $(time() - t1) s")
#         println(sens)
#         push!(t, sens)
#         save_name = savename(analysisDict)
#         CSV.write(scriptsdir("0nu/Bayes_hist_models/informative_priors/Cimrman_J40_ND_cuts/sensitivities_$(save_name)_$(rand(1:100000)).csv"), DataFrame(thalf= t))
#     catch
#         @warn "failed fit" 
#         continue
#     end
# end


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

data_hist = [get_pseudo_spectrum(b) for b in bkg_hist] |> sum 
# data_hist = sample_histogram(bkg_hist) 

# f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)


my_likelihood = make_hist_likelihood_uniform(
    data_hist,
    signal_hist_normed,
    bkg_hist_normed
)

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


burnin = MCMCMultiCycleBurnin(max_ncycles = 30, nsteps_final=1000)
mcmcalgo = MetropolisHastings(
    weighting = RepetitionWeighting(),
    tuning = AdaptiveMHTuning()
)


posterior = PosteriorMeasure(my_likelihood, prior)

@time logdensityof(my_likelihood, (As = 1e-5, Ab = fill(1e-5, length(bkg_hist))))

samples, _ = bat_sample(
    posterior, 
    MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = 5*10^4, nchains = 4),
)


binned_unshaped_samples, f_flatten = bat_transform(Vector, samples)
nDataPoints = integral(data_hist)
muS = [par[1] for par in binned_unshaped_samples.v]

@show exp_mu_signal_90 = quantile( muS,0.9) * nDataPoints 
Na = SNparams["Nₐ"]
m = SNparams["foilMass"] * SNparams["a"]
t = SNparams["tYear"]
W = SNparams["W"]
eff= lookup(signal, Bin_low, Bin_high-50)
Thalf = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90

Plots.plot(samples, size = (2000, 1600))
Plots.plot(samples, :As)
# savefig("notebooks/Sensitivity_exotic_Bayes_nDim/samples_As.png")

bin_width_new = 1#analysisDict[:bin_width]
xs = Bin_low:bin_width_new:Bin_high-bin_width_new/2
# b_amps = [x for x in mean(samples)][2:end] 
amps = mean(samples) 

# ys = [f2(amps, x) for x in xs]
# Plots.plot(data_hist, label = "data",fillrange = 1e-5, fillcolor = :match)
# Plots.plot!(xs, ys .* integral(data_hist) .* bin_width_new * bin_width, linewidth = 4,  label = "fit")
# # # Plots.plot!(xlims = (2000, 3200), ylims = (1e-5, 100))


# sig_amps = (As = quantile( muS,0.9), Ab = amps[2])
# ys_signal = [f2(sig_amps, x) for x in xs]
# Plots.plot!(xs, ys_signal .* integral(data_hist) .* bin_width_new * bin_width, linewidth = 4, label = "signal fit")
# Plots.plot!(yscale = :log10, ylims = (1.7e-2, 1e4))
# Plots.plot!(ylims = (0, 2))
# savefig("scripts/0nu/Bayes_hist_models/fit.svg")

# # SensitivityModule.stackedhist(bkg_hist)


using CairoMakie

colors = Makie.wong_colors()
f = Figure(size = (1200,800), fontsize = 26)
a = Axis(f[1,1], xlabel = L"\Theta_1", ylabel = L"Posterior density $$", title = L"Posterior density of signal fraction $\Theta_1$ \n for a single pseudo-experiment", limits = (0,nothing, 0, nothing))
CairoMakie.stephist!(a, muS, bins = range(0, maximum(muS)*1.1, length = 80), color = colors[1], label = "Posterior density", linewidth = 4, normalization = :pdf)
CairoMakie.vlines!(a, [quantile(muS, 0.9)], color = :red, linewidth = 4, linestyle= :dash, label = L"$\Theta_1$= %$(round(quantile(muS, 0.9), sigdigits = 3)) (90% C.I. upper limit)")
axislegend(a; position = :rt, patchsize = (45,35), labelsize = 20)
save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/nu0bb_radon1_posterior_signal_fraction.png"), f, px_per_unit = 4)
f


let 
    bkg_amps = amps[2] ./ sum(amps[2])
    labels = [L"2\nu\beta\beta", "Internal", "Radon", "Bi-210", "Detector", "External"]
    f = Figure(size = (1200,800), fontsize = 26)
    a = Axis(f[1,1], xlabel = L"$E_{sum}$ (keV)", ylabel = "Counts / 17.5kg.yr", title = "Fitted spectrum of a single pseudo-experiment", limits = (300,3500, 0, nothing) )
    h_total = sum(bkg_hist)
    n_total = integral(h_total)
    CairoMakie.stephist!(a, h_total, color = :black, label = "Total pseudo-background", linewidth = 4)
    
    h_fit = [ a*normalize(h, width = false)*n_total for (h, a) in zip(bkg_hist, bkg_amps)]

    for i in 1:length(h_fit[1:end])
        CairoMakie.stephist!(a, h_fit[i], color = colors[i], label = labels[i], linewidth = 4)
    end
    CairoMakie.hist!(a, sum(h_fit), label = "Fitted total", color = (colors[1], 0.4), )
    
    axislegend(a; position = :rt, patchsize = (45,35), labelsize = 20)
    f
end


begin 
    bin_width_new = 1
    xs = Bin_low:bin_width_new:Bin_high-bin_width_new/2
    # b_amps = [x for x in mean(samples)][2:end] 
    amps = mean(samples) 

    ys = [f2(amps, x) for x in xs]

    labels = [L"2\nu\beta\beta", "Internal", "Radon", "Bi-210", "Detector", "External"]
    f = Figure(size = (1200,800), fontsize = 26)
    a = Axis(f[1,1], xlabel = L"$E_{sum}$ (keV)", ylabel = L"Counts / 17.5kg.yr $$", title = "Fitted spectrum of a single pseudo-experiment", limits = (300,3500, 0, nothing) )
    h_total = sum(bkg_hist)
    n_total = integral(h_total)


    # CairoMakie.stephist!(a, h_total, color = :black, label = "Total pseudo-background", linewidth = 4)
    CairoMakie.scatter!(a, midpoints(binedges(h_total)), bincounts(h_total), label = "Pseudo-data", color = :black, markersize = 10)
    CairoMakie.errorbars!(a, h_total, color = :black, whiskerwidth=7, clamp_bincounts=true)

    CairoMakie.lines!(xs, ys .* integral(data_hist) .* bin_width_new * bin_width, linewidth = 2,  label = "total fit", linestyle= :dash)

    bkg_amps = amps[2] ./ sum(amps[2])
    h_fit = [ a*normalize(h, width = false)*n_total for (h, a) in zip(bkg_hist, bkg_amps)]
    for i in 1:length(h_fit[1:end])
        CairoMakie.stephist!(a, h_fit[i], color = colors[i], label = labels[i], linewidth = 4)
    end
    axislegend(a; position = :rt, patchsize = (45,35), labelsize = 20)

    xs_new = Bin_low:bin_width:Bin_high
    ys_new = [f2(amps, x) for x in xs_new]
    ys_total_fit = ys_new .* integral(data_hist) .* bin_width 
    residuals = bincounts(data_hist) .- ys_total_fit

    a2 = Axis(f[2,1], xlabel = L"$E_{sum}$ (keV)", ylabel = L"pseudo data - fit$$", limits = (300,3500, minimum(residuals)/2, maximum(residuals)*2  ) )
    CairoMakie.scatter!(a2, xs_new, residuals, color = :black, markersize = 10, label = "pseudo-data - fit")
    rowsize!(f.layout, 2, Relative(0.3))
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/nu0bb_radon1_fit.png"), f, px_per_unit = 4)
    f
end