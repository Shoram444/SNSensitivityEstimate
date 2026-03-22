using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using UnROOT
using BAT, SpecialFunctions, DensityInterface, IntervalSets, ValueShapes
using StatsPlots, ColorSchemes
# using LaTeXStrings
using MakieTeX

println("loaded pkgs")

include(datadir("sims/final_phd/data_sims_mul95_alpha_veto_curved_pcd/load_nd_sim_files.jl"))
include(scriptsdir("phd_final/08_data/bayes_fit_models.jl"))
include(scriptsdir("phd_final/08_data/helper_functions.jl"))

data_dir = datadir("data/final_phd/data_no_gamma_alpha_veto_nonan_fixed_us_window")
f1 = ROOTFile(joinpath(data_dir, "phase1_data.root"))
d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

f2 = ROOTFile(joinpath(data_dir, "phase2_data.root"))
d2 = LazyTree(f2, "tree", keys(f2["tree"])) |> DataFrame

f3 = ROOTFile(joinpath(data_dir, "phase3_data.root"))
d3 = LazyTree(f3, "tree", keys(f3["tree"])) |> DataFrame

var_data = :sumE
var_simu = :sumEsimu
binning = (500:100:3000)
fwhm = 0.14
data_cuts = Dict(
    :sumE => (700, 3000),
    :e1 => (350, 5000),
    :e2 => (350, 5000),
    :dy => (0, 100),
    :dz => (0, 100),
    :trackLength1 => (300, 1100),
    :trackLength2 => (300, 1100),
    :deltaCaloTime => (0, 1.5),
    :phi => (0,180)
)


simu_roi_1 = (
    sumEsimu = data_cuts[:sumE],
    simulatedEnergy1 = data_cuts[:e1],
    simulatedEnergy2 = data_cuts[:e2],
    dy = data_cuts[:dy],
    dz = data_cuts[:dz],
    trackLength1 = data_cuts[:trackLength1],
    trackLength2 = data_cuts[:trackLength2],
    # deltaCaloTime = data_cuts[:deltaCaloTime],
    # phi = data_cuts[:phi]
)
simu_roi_2 = (
    sumE = data_cuts[:sumE],
    reconstructedEnergy1 = data_cuts[:e1],
    reconstructedEnergy2 = data_cuts[:e2],
    dy = data_cuts[:dy],
    dz = data_cuts[:dz],
    trackLength1 = data_cuts[:trackLength1],
    trackLength2 = data_cuts[:trackLength2],
    # deltaCaloTime = data_cuts[:deltaCaloTime],
    # phi = data_cuts[:phi]
)

e1_min, e1_max = data_cuts[:e1]

data_roi_1 = (
    sumE = data_cuts[:sumE],
    reconstructedEnergy1 = data_cuts[:e1],
    reconstructedEnergy2 = data_cuts[:e2],
    dy = data_cuts[:dy],
    dz = data_cuts[:dz],
    trackLength1 = data_cuts[:trackLength1],
    trackLength2 = data_cuts[:trackLength2],
    deltaCaloTime = data_cuts[:deltaCaloTime],
    phi = data_cuts[:phi]
)

simu_varNames_1 = keys(simu_roi_1) |> collect .|> string
simu_p1_processes_roi_1 = load_sim_processes_ND(p1_dir, simu_varNames_1, simu_roi_1; timeMeas = p1_duration_seconds)
simu_p2_processes_roi_1 = load_sim_processes_ND(p2_dir, simu_varNames_1, simu_roi_1; timeMeas = p2_duration_seconds)
simu_p3_processes_roi_1 = load_sim_processes_ND(p3_dir, simu_varNames_1, simu_roi_1; timeMeas = p3_duration_seconds)

processes_old = load_data_processes( datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), "sumE"; roi = simu_roi_2 )
neutron_data_dir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/neutrons_jan_2026/")

include(joinpath(neutron_data_dir, "read_neutrons_1D.jl"))
neutron_processes = load_neutron_process1D(neutron_data_dir, "sumE"; roi=simu_roi_2)
for p in neutron_processes
    set_activity!(p, p.activity / 3)
    println("Neutron process: ", p.isotopeName, " nTotalSim: ", p.nTotalSim, " activity: ", p.activity)
end
backgrounds = [
    "Bi214_foil_bulk", 
    "Tl208_foil_bulk", 
    "Bi210_wire_surface",     
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
    "K40_hall_bulk"
]
background = [get_process(b, processes_old) |> first for b in backgrounds]
background = vcat(background, neutron_processes)


d1_roi_1 = filter_data!(d1, data_roi_1)
d2_roi_1 = filter_data!(d2, data_roi_1)
d3_roi_1 = filter_data!(d3, data_roi_1)
    

h1_data_p1 = Hist1D(getproperty(d1_roi_1, var_data); binedges = binning) #* inv(p1_duration_days)
h1_data_p2 = Hist1D(getproperty(d2_roi_1, var_data); binedges = binning) #* inv(p2_duration_days)
h1_data_p3 = Hist1D(getproperty(d3_roi_1, var_data); binedges = binning) #* inv(p3_duration_days)



p1_K40_process = get_process("K40_foil_bulk", simu_p1_processes_roi_1) |> first
p1_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p1_processes_roi_1) |> first
p1_bb_process = get_process("bb_foil_bulk", simu_p1_processes_roi_1) |> first
p1_radon_process = get_process("Bi214_wire_surface", simu_p1_processes_roi_1) |> first

p2_K40_process = get_process("K40_foil_bulk", simu_p2_processes_roi_1) |> first
p2_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p2_processes_roi_1) |> first
p2_bb_process = get_process("bb_foil_bulk", simu_p2_processes_roi_1) |> first
p2_radon_process = get_process("Bi214_wire_surface", simu_p2_processes_roi_1) |> first

p3_K40_process = get_process("K40_foil_bulk", simu_p3_processes_roi_1) |> first
p3_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p3_processes_roi_1) |> first
p3_bb_process = get_process("bb_foil_bulk", simu_p3_processes_roi_1) |> first
p3_radon_process = get_process("Bi214_wire_surface", simu_p3_processes_roi_1) |> first


simu_h_p1_K40 = get_simu_histo(p1_K40_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max)) 
simu_h_p1_Pa234m = get_simu_histo(p1_Pa234m_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))
simu_h_p1_bb = get_simu_histo(p1_bb_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))
simu_h_p1_radon = get_simu_histo(p1_radon_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))

simu_h_p2_K40 = get_simu_histo(p2_K40_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max)) 
simu_h_p2_Pa234m = get_simu_histo(p2_Pa234m_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))
simu_h_p2_bb = get_simu_histo(p2_bb_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))
simu_h_p2_radon = get_simu_histo(p2_radon_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))

simu_h_p3_K40 = get_simu_histo(p3_K40_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max)) 
simu_h_p3_Pa234m = get_simu_histo(p3_Pa234m_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))
simu_h_p3_bb = get_simu_histo(p3_bb_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))
simu_h_p3_radon = get_simu_histo(p3_radon_process, var_simu, fwhm; bins = binning, e_range = (e1_min, e1_max))


# h_flat_p2 = Hist1D(rand(Uniform(500, 3000),100000); binedges = binning)


function get_combined_bkg_hists(processes, binning)

    h_internal = Hist1D(;binedges = binning)
    h_radon = Hist1D(;binedges = binning)
    h_Bi210 = Hist1D(;binedges = binning)
    h_detector = Hist1D(;binedges = binning)
    h_external = Hist1D(;binedges = binning)

    internals = ["Bi214_foil_bulk", "Tl208_foil_bulk"]
    Bi210 = ["Bi210_wire_surface"]
    detectors = ["Bi214_PMT_glass_bulk","Bi214_PMT5_glass_bulk","K40_PMT_glass_bulk","K40_PMT5_glass_bulk","Tl208_PMT_glass_bulk","Tl208_PMT5_glass_bulk","Bi214_calo_8inch_scin_bulk","K40_calo_8inch_scin_bulk", ]

    for p in processes
        h = Hist1D(p.dataVector; binedges = binning)

        if p.isotopeName in internals && integral(h) >0
            h_internal = h_internal + h
        elseif p.isotopeName in Bi210 && integral(h) >0
            h_Bi210 = h_Bi210 + h
        elseif p.isotopeName in detectors && integral(h) >0
            h_detector = h_detector + h
        elseif integral(h) >0
            h_external = h_external + h
        end
    end
    return vcat(h_internal, h_Bi210, h_detector, h_external)
end

h_internal, h_Bi210, h_detector, h_external = get_combined_bkg_hists(background, binning)

proc_labels = [
    rich("2νββ"),
    rich(superscript("40"), "K"),
    rich(superscript("234m"), "Pa"),
    rich(superscript("214"), "Bi", " + ", superscript("208"), "Tl", " internal"),
    rich(superscript("214"), "Bi", "(radon)"),
    rich(superscript("210"), "Bi"),
    rich("detector"),
    rich("external"),
]

fit_results_free = []
fit_component_stds = Vector{Vector{Float64}}()

samples_by_phase = []

for p in 1:3
    println("Fitting phase $p")
    phase = p

    if phase == 1
        h_K40, h_Pa234m, h_bb, h_radon = simu_h_p1_K40, simu_h_p1_Pa234m, simu_h_p1_bb, simu_h_p1_radon
    elseif phase == 2
        h_K40, h_Pa234m, h_bb, h_radon = simu_h_p2_K40, simu_h_p2_Pa234m, simu_h_p2_bb, simu_h_p2_radon
    elseif phase == 3
        h_K40, h_Pa234m, h_bb, h_radon = simu_h_p3_K40, simu_h_p3_Pa234m, simu_h_p3_bb, simu_h_p3_radon
    end 


    bkg_hist = [
        h_bb,
        h_K40,
        h_Pa234m,
        h_internal,
        h_radon,
        h_Bi210,
        h_detector,
        h_external,
    ]

    data_hist = if phase == 1
        h1_data_p1
    elseif phase == 2
        h1_data_p2
    elseif phase == 3
        h1_data_p3
    end

    bkg_hist_normed = normalize.(bkg_hist, width = false)


    nData = sum(bincounts(data_hist))

    μ_K40_prior      = Uniform(0.0, 2*nData)
    μ_Pa_prior       = Uniform(0.0, 2*nData)
    μ_bb_prior       = Uniform(0.0, 2*nData)
    μ_radon_prior    = Uniform(0.0, 2*nData)
    μ_internal_prior = Uniform(0.0, 2*nData)
    μ_Bi210_prior    = Uniform(0.0, 2*nData)
    μ_detector_prior = Uniform(0.0, 2*nData)
    μ_external_prior = Uniform(0.0, 2*nData)

    prior_free = NamedTupleDist(
        μ_bb = μ_bb_prior,
        μ_K40 = μ_K40_prior,
        μ_Pa = μ_Pa_prior,
        μ_internal = μ_internal_prior,
        μ_radon = μ_radon_prior,
        μ_Bi210 = μ_Bi210_prior,
        μ_detector = μ_detector_prior,
        μ_external = μ_external_prior
    )
    


    likelihood_free = make_hist_likelihood_free(data_hist, bkg_hist_normed)
    posterior_free = PosteriorMeasure(likelihood_free, prior_free)

    burnin = MCMCMultiCycleBurnin(max_ncycles=30, nsteps_final=3000)
    mcmcalgo = RandomWalk()

    samples, _ = bat_sample(
        posterior_free,
        TransformedMCMC(
            proposal = mcmcalgo,
            burnin=burnin,
            nsteps=5*10^5,
            nchains=4
        )
    ) 

    push!(samples_by_phase, samples)

    Plots.plot(samples, vsel = (1:10),size = (8000, 8000), xrotation=45, thickness_scaling = 1.5,fontsize = 90, margin= 10Plots.mm, dpi = 120)
    savefig(current(), scriptsdir("phd_final/08_data/figs/unconstrained_parameter_cor_phase$(phase).png"))

    post_mean = mean(samples)

    μ_K40 = post_mean.μ_K40
    μ_Pa = post_mean.μ_Pa
    μ_bb = post_mean.μ_bb
    μ_radon = post_mean.μ_radon
    μ_internal = post_mean.μ_internal
    μ_Bi210 = post_mean.μ_Bi210
    μ_detector = post_mean.μ_detector
    μ_external = post_mean.μ_external

    std_post = std(samples)
    σ_bb = std_post.μ_bb
    σ_K40 = std_post.μ_K40
    σ_Pa = std_post.μ_Pa
    σ_internal = std_post.μ_internal
    σ_radon = std_post.μ_radon
    σ_Bi210 = std_post.μ_Bi210
    σ_detector = std_post.μ_detector
    σ_external = std_post.μ_external

    μ_vec = [
        μ_bb,
        μ_K40,
        μ_Pa,
        μ_internal,
        μ_radon,
        μ_Bi210,
        μ_detector,
        μ_external
    ]

    σ_vec = [
        σ_bb,
        σ_K40,
        σ_Pa,
        σ_internal,
        σ_radon,
        σ_Bi210,
        σ_detector,
        σ_external,
    ]

    function build_fit_histograms(bkg_hists_normed, μ)
        hists = Hist1D[]
        for i in 1:length(bkg_hists_normed)
            push!(hists, μ[i] * bkg_hists_normed[i])
        end
        return hists
    end

    fit_hists = build_fit_histograms(bkg_hist_normed, μ_vec)
    total_fit = sum(fit_hists)

    push!(
        fit_results_free, 
        fit_result(
            phase,
            data_hist,
            fit_hists,
            post_mean,
            proc_labels
        )
    )
    push!(fit_component_stds, σ_vec)


    activities = [
        μ_bb,
        μ_K40,
        μ_Pa,
        μ_internal,
        μ_radon,
        μ_Bi210,
        μ_detector,
        μ_external,
    ]


    plot_fit(
        phase = phase,
        data_hist = data_hist,
        fitted_hists = fit_hists,
        fit_params = post_mean,
        component_count_uncertainties = σ_vec,
        show_component_uncertainties = true,
        outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_phase$(phase)_free.png"),
        component_labels = proc_labels,
        main_limits = (0, 4000, 0, nothing),
        ratio_limits = (0, 4000, 0.2, 1.8),
        blinded_roi = (2700, 3000),
        chi2_nparams = length(μ_vec),
        include_component_counts = true,
    )

end


################
## COMBINE PHASES
################

full_data = sum([r.data_hist for r in fit_results_free])


fit_p1 = fit_results_free[1].fitted_hists
fit_p2 = fit_results_free[2].fitted_hists
fit_p3 = fit_results_free[3].fitted_hists

fit_combined = sum([fit_p1, fit_p2, fit_p3])
fit_total = sum([sum(fit_p1), sum(fit_p2), sum(fit_p3)])
combined_component_stds = sqrt.(fit_component_stds[1].^2 .+ fit_component_stds[2].^2 .+ fit_component_stds[3].^2)



f_full = plot_fit(
    phase = 0,
    data_hist = full_data,
    fitted_hists = fit_combined,
    fit_params = nothing,
    component_count_uncertainties = combined_component_stds,
    show_component_uncertainties = true,
    outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_combined_free_lines.png"),
    component_labels = proc_labels,
    main_limits = (0, 4000, 0, nothing),
    ratio_limits = (0, 4000, 0.2, 1.8),
    blinded_roi = (2700, 3000),
    chi2_nparams = length(fit_combined),
    title = "Combined Fit phase 1+2+3",
    include_component_counts = true,
    stacked = false,
    filled = false,
)





#### plot samples corrplot


p0 = StatsPlots.plot(
    samples_by_phase[1],
    margin = 2Plots.mm,
    size = (3600, 3600),
    dpi =60,
    thickness_scaling = 2,
)

p1 = StatsPlots.plot(
    samples_by_phase[1],
    vsel = (:μ_bb, :μ_radon),
    margin = 10Plots.mm,
    figure_title = "Posterior samples of μ_bb vs μ_radon for phase 1",
)
savefig(p1, scriptsdir("phd_final/08_data/figs/unconstrained_parameter_cor_phase1_bb_radon.png"))

p2 = StatsPlots.plot(
    samples_by_phase[1],
    vsel = (:μ_radon, :μ_internal),
    margin = 10Plots.mm,
    figure_title = "Posterior samples of μ_radon vs μ_internal for phase 1",
)
savefig(p2, scriptsdir("phd_final/08_data/figs/unconstrained_parameter_cor_phase1_radon_internal.png"))