using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using UnROOT
using BAT, SpecialFunctions, DensityInterface, IntervalSets, ValueShapes
using StatsPlots, ColorSchemes, Measurements
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
    deltaCaloTime = data_cuts[:deltaCaloTime],
    phi = data_cuts[:phi]
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

processes_old = load_data_processes( datadir("sims/final_phd/data_sims_mul95_alpha_veto_curved_pcd/old_sims"), "sumE"; roi = simu_roi_2 )
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


function get_process_info(processes, isotopeName; var_simu = :sumEsimu, fwhm = 0.14, bins = binning, e_range = (e1_min, e1_max), duration = p1_duration_seconds)
    proc = get_process(isotopeName, processes) |> first
    h = get_simu_histo(proc, var_simu, fwhm; bins = bins, e_range = e_range) 
    efficiency = integral(h) / proc.nTotalSim
    activity = proc.activity
    amount = proc.amount
    n_exp = activity * duration * efficiency * amount
    return h, efficiency, activity, n_exp, proc
end

simu_h1_p1_K40, eff_p1_K40, a_p1_K40, n_exp_p1_K40, p1_K40_process = get_process_info(simu_p1_processes_roi_1, "K40_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p1_duration_seconds)
simu_h1_p1_Pa234m, eff_p1_Pa234m, a_p1_Pa234m, n_exp_p1_Pa234m, p1_Pa234m_process = get_process_info(simu_p1_processes_roi_1, "Pa234m_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p1_duration_seconds)
simu_h1_p1_bb, eff_p1_bb, a_p1_bb, n_exp_p1_bb, p1_bb_process = get_process_info(simu_p1_processes_roi_1, "bb_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p1_duration_seconds)
simu_h1_p1_radon, eff_p1_radon, a_p1_radon, n_exp_p1_radon, p1_radon_process = get_process_info(simu_p1_processes_roi_1, "Bi214_wire_surface"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p1_duration_seconds)

simu_h1_p2_K40, eff_p2_K40, a_p2_K40, n_exp_p2_K40, p2_K40_process = get_process_info(simu_p2_processes_roi_1, "K40_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p2_duration_seconds)
simu_h1_p2_Pa234m, eff_p2_Pa234m, a_p2_Pa234m, n_exp_p2_Pa234m, p2_Pa234m_process = get_process_info(simu_p2_processes_roi_1, "Pa234m_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p2_duration_seconds)
simu_h1_p2_bb, eff_p2_bb, a_p2_bb, n_exp_p2_bb, p2_bb_process = get_process_info(simu_p2_processes_roi_1, "bb_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p2_duration_seconds)
simu_h1_p2_radon, eff_p2_radon, a_p2_radon, n_exp_p2_radon, p2_radon_process = get_process_info(simu_p2_processes_roi_1, "Bi214_wire_surface"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p2_duration_seconds)

simu_h1_p3_K40, eff_p3_K40, a_p3_K40, n_exp_p3_K40, p3_K40_process = get_process_info(simu_p3_processes_roi_1, "K40_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p3_duration_seconds)
simu_h1_p3_Pa234m, eff_p3_Pa234m, a_p3_Pa234m, n_exp_p3_Pa234m, p3_Pa234m_process = get_process_info(simu_p3_processes_roi_1, "Pa234m_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p3_duration_seconds)
simu_h1_p3_bb, eff_p3_bb, a_p3_bb, n_exp_p3_bb, p3_bb_process = get_process_info(simu_p3_processes_roi_1, "bb_foil_bulk"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p3_duration_seconds)
simu_h1_p3_radon, eff_p3_radon, a_p3_radon, n_exp_p3_radon, p3_radon_process = get_process_info(simu_p3_processes_roi_1, "Bi214_wire_surface"; var_simu = var_simu, fwhm = fwhm, bins = binning, e_range = (e1_min, e1_max), duration = p3_duration_seconds)


function get_combined_bkg_hists(processes, binning)

    h_internal = Hist1D(;binedges = binning)
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

##### Constraining parameters we know of

function make_hist_likelihood_constrained(
    data_hist,
    bkg_hists;
    C_K40,
    C_Pa,
    C_bb,
    C_radon,
    a_bb_fixed,
    # eff_bb,
    a_radon_fixed,
    # eff_radon,
)

    data_counts = bincounts(data_hist)
    bkg_counts = [bincounts(h) for h in bkg_hists]

    nbins = length(data_counts)

    logfuncdensity(function(p)
        μ_bb = C_bb * a_bb_fixed * p.eff_bb
        μ_K40 = C_K40 * p.a_K40 * p.eff_K40
        μ_Pa = C_Pa * p.a_Pa * p.eff_Pa
        μ_internal = p.μ_internal
        μ_radon = C_radon * a_radon_fixed * p.eff_radon
        μ_Bi210 = p.μ_Bi210
        μ_detector = p.μ_detector
        μ_external = p.μ_external

        μ_vec = (μ_bb, μ_K40, μ_Pa, μ_internal, μ_radon, μ_Bi210, μ_detector, μ_external)

        ll = 0.0

        for k in 1:nbins
            λk = 0.0

            for i in 1:length(bkg_hists)
                λk += μ_vec[i] * bkg_counts[i][k]
            end

            λk = max(λk, eps())
            ll += log_pdf_poisson(λk, data_counts[k])
        end

        return ll
    end)
end

function lognormal_from_mean_relerr(mean, relerr)
    σ_log = sqrt(log(1 + relerr^2))
    μ_log = log(mean) - 0.5 * σ_log^2
    return LogNormal(μ_log, σ_log)
end

fit_results_informed = []
fit_component_stds_informed = Vector{Vector{Float64}}()
radon_activity_by_phase = Dict(1 => 110.1e-3, 2 => 67.7e-3, 3 => 3.36e-3)
samples_by_phase = []

for p in 1:3
    println("Fitting informed phase $p")
    phase = p

    if phase == 1
        h_K40, h_Pa234m, h_bb, h_radon = simu_h1_p1_K40, simu_h1_p1_Pa234m, simu_h1_p1_bb, simu_h1_p1_radon
        eff_K40_ref, eff_Pa_ref, eff_bb_ref, eff_radon_ref = eff_p1_K40, eff_p1_Pa234m, eff_p1_bb, eff_p1_radon
        duration = p1_duration_seconds
        p_K40_proc, p_Pa_proc, p_bb_proc, p_radon_proc = p1_K40_process, p1_Pa234m_process, p1_bb_process, p1_radon_process
    elseif phase == 2
        h_K40, h_Pa234m, h_bb, h_radon = simu_h1_p2_K40, simu_h1_p2_Pa234m, simu_h1_p2_bb, simu_h1_p2_radon
        eff_K40_ref, eff_Pa_ref, eff_bb_ref, eff_radon_ref = eff_p2_K40, eff_p2_Pa234m, eff_p2_bb, eff_p2_radon
        duration = p2_duration_seconds
        p_K40_proc, p_Pa_proc, p_bb_proc, p_radon_proc = p2_K40_process, p2_Pa234m_process, p2_bb_process, p2_radon_process
    else
        h_K40, h_Pa234m, h_bb, h_radon = simu_h1_p3_K40, simu_h1_p3_Pa234m, simu_h1_p3_bb, simu_h1_p3_radon
        eff_K40_ref, eff_Pa_ref, eff_bb_ref, eff_radon_ref = eff_p3_K40, eff_p3_Pa234m, eff_p3_bb, eff_p3_radon
        duration = p3_duration_seconds
        p_K40_proc, p_Pa_proc, p_bb_proc, p_radon_proc = p3_K40_process, p3_Pa234m_process, p3_bb_process, p3_radon_process
    end

    data_hist = if phase == 1
        h1_data_p1
    elseif phase == 2
        h1_data_p2
    else
        h1_data_p3
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
    bkg_hist_normed = normalize.(bkg_hist, width = false)
    nData = sum(bincounts(data_hist))

    C_K40 = duration * p_K40_proc.amount
    C_Pa = duration * p_Pa_proc.amount
    C_bb = duration * p_bb_proc.amount
    C_radon = duration * p_radon_proc.amount

    a_bb_fixed = p_bb_proc.activity
    # eff_bb_fixed = eff_bb_ref
    a_radon_fixed = radon_activity_by_phase[phase]
    # eff_radon_fixed = eff_radon_ref

    prior_informed = NamedTupleDist(
        eff_bb = lognormal_from_mean_relerr(eff_bb_ref, 0.1),
        eff_radon = lognormal_from_mean_relerr(eff_radon_ref, 0.1),
        a_K40 = Uniform(0.0, 2 * nData),
        eff_K40 = lognormal_from_mean_relerr(eff_K40_ref, 0.3),
        a_Pa = Uniform(0.0, 2 * nData),
        eff_Pa = lognormal_from_mean_relerr(eff_Pa_ref, 0.3),
        μ_internal = Uniform(0.0, 2 * nData),
        μ_Bi210 = Uniform(0.0, 2 * nData),
        μ_detector = Uniform(0.0, 2 * nData),
        μ_external = Uniform(0.0, 2 * nData),
    )

    likelihood_informed = make_hist_likelihood_constrained(
        data_hist,
        bkg_hist_normed;
        C_K40 = C_K40,
        C_Pa = C_Pa,
        C_bb = C_bb,
        C_radon = C_radon,
        a_bb_fixed = a_bb_fixed,
        # eff_bb_fixed = eff_bb_fixed,
        a_radon_fixed = a_radon_fixed,
        # eff_radon_fixed = eff_radon_fixed,
    )
    posterior_informed = PosteriorMeasure(likelihood_informed, prior_informed)

    burnin = MCMCMultiCycleBurnin(max_ncycles = 30, nsteps_final = 3000)
    mcmcalgo = RandomWalk()

    samples, _ = bat_sample(
        posterior_informed,
        TransformedMCMC(
            proposal = mcmcalgo,
            burnin = burnin,
            nsteps = 5 * 10^5,
            nchains = 4,
        )
    )

    # StatsPlots.plot(samples, vsel = (1:10),size = (8000, 8000), xrotation=45,thickness_scaling = 1.6, margin= 1Plots.mm, dpi = 100)
    # savefig(current(), scriptsdir("phd_final/08_data/figs/informed_parameter_cor_phase$(phase).png"))

    push!(samples_by_phase, samples)

    post_mean = mode(samples) #mean(samples)
    std_post = std(samples)

    mm = map(measurement, post_mean, std_post)

    μ_bb = C_bb * a_bb_fixed * mm.eff_bb
    μ_K40 = C_K40 * mm.a_K40 * mm.eff_K40
    μ_Pa = C_Pa * mm.a_Pa * mm.eff_Pa
    μ_internal = mm.μ_internal
    μ_radon = C_radon * a_radon_fixed * mm.eff_radon
    μ_Bi210 = mm.μ_Bi210
    μ_detector = mm.μ_detector
    μ_external = mm.μ_external

    σ_bb = Measurements.uncertainty(μ_bb)
    σ_K40 = Measurements.uncertainty(μ_K40)
    σ_Pa = Measurements.uncertainty(μ_Pa)
    σ_internal = Measurements.uncertainty(μ_internal)
    σ_radon = Measurements.uncertainty(μ_radon)
    σ_Bi210 = Measurements.uncertainty(μ_Bi210)
    σ_detector = Measurements.uncertainty(μ_detector)
    σ_external = Measurements.uncertainty(μ_external)

    μ_vec = [
        μ_bb,
        μ_K40,
        μ_Pa,
        μ_internal,
        μ_radon,
        μ_Bi210,
        μ_detector,
        μ_external,
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

    n_vec = measurement.(μ_vec, σ_vec)

    fit_hists = build_fit_histograms(bkg_hist_normed, Measurements.value.(μ_vec))

    push!(
        fit_results_informed,
        fit_result(
            phase,
            data_hist,
            fit_hists,
            post_mean,
            n_vec,
            proc_labels,
        )
    )
    push!(fit_component_stds_informed, σ_vec)

    plot_fit(
        phase = phase,
        data_hist = data_hist,
        fitted_hists = fit_hists,
        fit_params = post_mean,
        component_count_uncertainties = n_vec,
        show_component_uncertainties = true,
        outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_phase$(phase)_informed.png"),
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

full_data = sum([r.data_hist for r in fit_results_informed])

fit_p1 = fit_results_informed[1].fitted_hists
fit_p2 = fit_results_informed[2].fitted_hists
fit_p3 = fit_results_informed[3].fitted_hists

fit_combined = sum([fit_p1, fit_p2, fit_p3])
combined_component_stds_informed = sqrt.(fit_component_stds_informed[1].^2 .+ fit_component_stds_informed[2].^2 .+ fit_component_stds_informed[3].^2)

f_full_informed = plot_fit(
    phase = 0,
    data_hist = full_data,
    fitted_hists = fit_combined,
    fit_params = nothing,
    component_count_uncertainties = combined_component_stds_informed,
    show_component_uncertainties = true,
    outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_combined_informed.png"),
    component_labels = proc_labels,
    main_limits = (0, 4000, 0, nothing),
    ratio_limits = (0, 4000, 0.2, 1.8),
    blinded_roi = (2700, 3000),
    chi2_nparams = length(fit_combined),
    title = "Combined informed fit phase 1+2+3",
    include_component_counts = true,
    stacked = false,
    filled = false,
)


# Plot the priors and posteriors of each component for each phase
fig_bb_prior, _, _ = plot_lognormal_prior_comparison(;
    mean = eff_p1_bb,
    relerr = 0.1,
    title = nothing,
    lognormal_title = "LogNormal prior for 2νββ efficiency",
    lognormal_xlabel = "ε",
    ylabel = "Density",
    figure_size = (1200, 350),
    fontsize = 28,
    n_sigma = 4,
    n_points = 1200,
    px_per_unit = 3,
    outputpath = scriptsdir("phd_final/08_data/figs/prior_eff_bb.png"),
)

fig_radon_prior, _, _ = plot_lognormal_prior_comparison(;
    mean = eff_p1_radon,
    relerr = 0.1,
    title = nothing,
    lognormal_title = "LogNormal prior for radon efficiency",
    lognormal_xlabel = "ε",
    ylabel = "Density",
    figure_size = (1200, 350),
    fontsize = 28,
    n_sigma = 4,
    n_points = 1200,
    px_per_unit = 3,
    outputpath = scriptsdir("phd_final/08_data/figs/prior_eff_radon.png"),
)

fig_K40_prior, _, _ = plot_lognormal_prior_comparison(;
    mean = eff_p1_K40,
    relerr = 0.3,
    title = nothing,
    lognormal_title = "LogNormal prior for K40 efficiency",
    lognormal_xlabel = "ε",
    ylabel = "Density",
    figure_size = (1200, 350),
    fontsize = 28,
    n_sigma = 4,
    n_points = 1200,
    px_per_unit = 3,
    outputpath = scriptsdir("phd_final/08_data/figs/prior_eff_K40.png"),
)

fig_Pa234m_prior, _, _ = plot_lognormal_prior_comparison(;
    mean = eff_p1_Pa234m,
    relerr = 0.3,
    title = nothing,
    lognormal_title = "LogNormal prior for Pa234m efficiency",
    lognormal_xlabel = "ε",
    ylabel = "Density",
    figure_size = (1200, 350),
    fontsize = 28,
    n_sigma = 4,
    n_points = 1200,
    px_per_unit = 3,
    outputpath = scriptsdir("phd_final/08_data/figs/prior_eff_Pa234m.png"),
)



let 
    f = Figure(size = (1200, 350), fontsize = 28, figure_padding = 35)
    nData = Int(h1_data_p1 |> integral)

    a = Axis(f[1,1], title = "Uniform prior for μ", xlabel = "μ", ylabel = "Density", limits = (0, 2*nData, 0, inv(2*nData)*1.5))
    
    x = range(0, stop = 2*nData, length = 1000) 
    y = pdf.(Uniform(0, 2*nData), x)
    band!(a, x, zero.(y), y, color = ("#003865", 0.4))
    lines!(a, x, y, color ="#003865", label = "Uniform(0, 2*nData)",linewidth = 4)
    axislegend(a)
    save(scriptsdir("phd_final/08_data/figs/prior_mu_others.png"), f)
    f
end


function get_activity_from_n(n, t, m, e)
    # n: expected counts
    # t: duration in seconds
    # m: mass in kg
    # e: efficiency (0-1)
    return n / (t * m * e)
end

n1_radon = measurement(256.7,22.2)
n2_radon = measurement(642.6,58.8)
n3_radon = measurement(16.9, 1.7)

a1_radon = get_activity_from_n(n1_radon, p1_duration_seconds, p1_radon_process.amount, fit_results_informed[1].params.eff_radon)
a2_radon = get_activity_from_n(n2_radon, p2_duration_seconds, p2_radon_process.amount, fit_results_informed[2].params.eff_radon)
a3_radon = get_activity_from_n(n3_radon, p3_duration_seconds, p3_radon_process.amount, fit_results_informed[3].params.eff_radon)
a_radon_combined = ( a1_radon * p1_duration_seconds + a2_radon * p2_duration_seconds + a3_radon * p3_duration_seconds ) / (p1_duration_seconds + p2_duration_seconds + p3_duration_seconds)

total_duration = p1_duration_seconds + p2_duration_seconds + p3_duration_seconds
mean_eff_pa = mean([fit_results_informed[1].params.eff_Pa, fit_results_informed[2].params.eff_Pa, fit_results_informed[3].params.eff_Pa])
n_pa = measurement(125.8, 72.7)
a_pa234m = get_activity_from_n(n_pa, total_duration, p1_Pa234m_process.amount, mean_eff_pa)

mean_eff_k40 = mean([fit_results_informed[1].params.eff_K40, fit_results_informed[2].params.eff_K40, fit_results_informed[3].params.eff_K40])
n_K40 = measurement(47.2,27.3)
a_K40 = get_activity_from_n(n_K40, total_duration, p1_K40_process.amount, mean_eff_k40)

p_Bi210_process = get_process("Bi210_wire_surface", background) |> first
n_bi210 = measurement(40.3, 16.6)
n_bi210_exp = total_duration * p_Bi210_process.amount * (h_Bi210 |> integral) / p_Bi210_process.nTotalSim * p_Bi210_process.activity
r_bi210 = n_bi210 / n_bi210_exp
a_bi210 = r_bi210 * p_Bi210_process.activity



#### 

with_theme(theme_latexfonts()) do 
    f = Figure(size = (1200, 700), fontsize = 28, figure_padding = 35)
    a = Axis(f[1,1], title = "2-electron spectrum of runs 2011 - 3743", xlabel = L"E_{sum} (\text{keV})", ylabel = L"Counts / $100\, \text{keV}$", limits = (0, 3500, 0, nothing))

    h_tot = h1_data_p1 + h1_data_p2 + h1_data_p3
    n1, n2, n3, ntot = Int.(integral.([h1_data_p1, h1_data_p2, h1_data_p3, h_tot]))
    CairoMakie.scatter!(a, h1_data_p1, color = colors[1], label = "Phase 1, n = $(n1)", markersize = 8)
    CairoMakie.errorbars!(a, h1_data_p1, color = colors[1], whiskerwidth = 10, label = "Phase 1, n = $(n1)")
    CairoMakie.scatter!(a, h1_data_p2, color = colors[2], label = "Phase 2, n = $(n2)", markersize = 8)
    CairoMakie.errorbars!(a, h1_data_p2, color = colors[2], whiskerwidth = 10, label = "Phase 2, n = $(n2)")
    CairoMakie.scatter!(a, h1_data_p3, color = colors[3], label = "Phase 3, n = $(n3)", markersize = 8)
    CairoMakie.errorbars!(a, h1_data_p3, color = colors[3], whiskerwidth = 10, label = "Phase 3, n = $(n3)")

    CairoMakie.scatter!(a, h_tot, color = :black, label = "Combined, n = $(ntot)", markersize = 8)
    CairoMakie.errorbars!(a, h_tot, color = :black, whiskerwidth = 10, label = "Combined, n = $(ntot)")

    CairoMakie.vspan!(a, 2700, 3000, color = :black, label = "Blinded ROI")

    axislegend(a, position = :rt, merge = true)
    save(scriptsdir("phd_final/08_data/figs/data_spectrum_phases.png"), f)
    f
end



p1_oms = CSV.File("data/data/final_phd/simrc/dead_om_list_phase1.csv") |>DataFrame
p2_oms = CSV.File("data/data/final_phd/simrc/dead_om_list_phase2.csv") |>DataFrame
p3_oms = CSV.File("data/data/final_phd/simrc/dead_om_list_phase3.csv") |>DataFrame

select!(p1_oms, [:om_num => :phase1_oms])
select!(p2_oms, [:om_num => :phase2_oms])
select!(p3_oms, [:om_num => :phase3_oms])

function _pad_with_missing(col, n)
    out = Vector{Union{Missing, eltype(col)}}(missing, n)
    out[1:length(col)] = col
    return out
end

max_len_oms = maximum((nrow(p1_oms), nrow(p2_oms), nrow(p3_oms)))

p_oms = DataFrame(
    phase1_oms = _pad_with_missing(p1_oms.phase1_oms, max_len_oms),
    phase2_oms = _pad_with_missing(p2_oms.phase2_oms, max_len_oms),
    phase3_oms = _pad_with_missing(p3_oms.phase3_oms, max_len_oms),
)


open(joinpath(scriptsdir("phd_final/08_data/"), "dead_om_list_combined.tex"), "w") do io
    pretty_table(
        io,
        p_oms,
        header = ["Phase 1 OM number", "Phase 2 OM number", "Phase 3 OM number"],
        backend = Val(:latex)
    )
end

sum(p3_oms.phase3_oms .< 519)





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
    vsel = (:eff_bb, :eff_radon),
    margin = 10Plots.mm,
    figure_title = "Posterior samples of ε_bb vs ε_radon for phase 1",
)
savefig(p1, scriptsdir("phd_final/08_data/figs/constrained_parameter_cor_phase1_bb_radon.png"))

p2 = StatsPlots.plot(
    samples_by_phase[1],
    vsel = (:eff_radon, :μ_internal),
    margin = 10Plots.mm,
    figure_title = "Posterior samples of ε_radon vs μ_internal for phase 1",
)
savefig(p2, scriptsdir("phd_final/08_data/figs/constrained_parameter_cor_phase1_radon_internal.png"))

