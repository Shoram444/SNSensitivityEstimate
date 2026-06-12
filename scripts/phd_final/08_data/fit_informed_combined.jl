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

begin
var_data = :sumE
var_simu = :sumEsimu
binning = (500:50:2700)
fwhm = 0.15
data_cuts = Dict(
    :sumE => (500, 2700),
    :e1 => (250, 5000),
    :e2 => (250, 5000),
    :dy => (0, 80),
    :dz => (0, 40),
    :trackLength1 => (300, 1100),
    :trackLength2 => (300, 1100),
    :deltaCaloTime => (0, 2.2),
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
    phi = data_cuts[:phi]
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

simu_var2 = var_simu == :sumEsimu ? :sumE : var_simu
processes_old = load_data_processes( datadir("sims/final_phd/data_sims_mul95_alpha_veto_curved_pcd/old_sims"), string(simu_var2); roi = simu_roi_2 )
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
# background = vcat(background, neutron_processes)


d1_roi_1 = filter_data!(d1, data_roi_1)
d2_roi_1 = filter_data!(d2, data_roi_1)
d3_roi_1 = filter_data!(d3, data_roi_1)
    

h1_data_p1 = Hist1D(getproperty(d1_roi_1, var_data); binedges = binning) #* inv(p1_duration_days)
h1_data_p2 = Hist1D(getproperty(d2_roi_1, var_data); binedges = binning) #* inv(p2_duration_days)
h1_data_p3 = Hist1D(getproperty(d3_roi_1, var_data); binedges = binning) #* inv(p3_duration_days)

# Combined data histogram
h1_data_combined = h1_data_p1 + h1_data_p2 + h1_data_p3

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
    rich(superscript("214"), "Bi", " (radon)"),
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

println("Fitting informed combined fit (all phases together)")

# Combined histograms
h_K40_combined = simu_h1_p1_K40 + simu_h1_p2_K40 + simu_h1_p3_K40
h_Pa234m_combined = simu_h1_p1_Pa234m + simu_h1_p2_Pa234m + simu_h1_p3_Pa234m
h_bb_combined = simu_h1_p1_bb + simu_h1_p2_bb + simu_h1_p3_bb
h_radon_combined = simu_h1_p1_radon + simu_h1_p2_radon + simu_h1_p3_radon

# Mean efficiencies
eff_K40_ref = mean([eff_p1_K40, eff_p2_K40, eff_p3_K40])
eff_Pa_ref = mean([eff_p1_Pa234m, eff_p2_Pa234m, eff_p3_Pa234m])
eff_bb_ref = mean([eff_p1_bb, eff_p2_bb, eff_p3_bb])
eff_radon_ref = mean([eff_p1_radon, eff_p2_radon, eff_p3_radon])

# Total duration
total_duration = p1_duration_seconds + p2_duration_seconds + p3_duration_seconds

# Combined processes (use phase 1 as reference but account for all phases in amounts)
K40_process = p1_K40_process
Pa234m_process = p1_Pa234m_process
bb_process = p1_bb_process
radon_process = p1_radon_process

# Use combined amounts
K40_process_combined = deepcopy(K40_process)
K40_process_combined.timeMeas = p1_K40_process.timeMeas + p2_K40_process.timeMeas + p3_K40_process.timeMeas

Pa234m_process_combined = deepcopy(Pa234m_process)
Pa234m_process_combined.timeMeas = p1_Pa234m_process.timeMeas + p2_Pa234m_process.timeMeas + p3_Pa234m_process.timeMeas

bb_process_combined = deepcopy(bb_process)
bb_process_combined.timeMeas = p1_bb_process.timeMeas + p2_bb_process.timeMeas + p3_bb_process.timeMeas

radon_process_combined = deepcopy(radon_process)
radon_process_combined.timeMeas = p1_radon_process.timeMeas + p2_radon_process.timeMeas + p3_radon_process.timeMeas

bkg_hist = [
    h_bb_combined,
    h_K40_combined,
    h_Pa234m_combined,
    h_internal,
    h_radon_combined,
    h_Bi210,
    h_detector,
    h_external,
]
bkg_hist_normed = normalize.(bkg_hist, width = false)
nData = sum(bincounts(h1_data_combined))

C_K40 = total_duration * K40_process_combined.amount
C_Pa = total_duration * Pa234m_process_combined.amount
C_bb = total_duration * bb_process_combined.amount
C_radon = total_duration * radon_process_combined.amount

a_bb_fixed = bb_process_combined.activity
a_radon_fixed = 53.1e-3  # 53.1 mBq/m converted to Bq/m

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
    h1_data_combined,
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

post_mean = mean(samples)
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

fit_hists = build_fit_histograms(bkg_hist_normed, Measurements.value.(μ_vec))

plot_style_colors = [
    "#0F4C81",  # 2νββ
    "#D17A22",  # 40K
    "#2F4858",  # 234mPa
    "#A23E2A",  # internal
    "#7A2E8C",  # radon
    "#2A9D8F",  # Bi210
    "#4A4E69",  # detector
    "#6C757D",  # external
]

fit_outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_informed_combined.png")
fit_pull_outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_informed_combined_pull.png")



p1 = plot_fit(
    phase = 0,
    data_hist = h1_data_combined,
    fitted_hists = fit_hists,
    fit_params = post_mean,
    component_count_uncertainties = μ_vec,
    show_component_uncertainties = true,
    title = "Combined Informed Fit all phases",
    outputpath = fit_outputpath,
    component_labels = proc_labels,
    main_limits = (0, 4000, 0, nothing),
    ratio_limits = (0, 4000, 0., 2.0),
    blinded_roi = (2700, 3000),
    chi2_nparams = length(μ_vec),
    include_component_counts = true,
    figure_size = (3400*0.75, 1350*0.75),
    chi2_text_pos = (0.79, 0.72),
    fontsize = 44,
    colors = plot_style_colors,
    component_linewidth = 5,
    data_markersize = 16,
    data_whiskerwidth = 12,
    total_fit_color = :red,
    total_fit_linewidth = 7,
    legend_tellheight = true,
    legend_width = 950,
    legend_nbanks = 2,
    legend_rowgap = 40
)


p2 = plot_fit_pull(
    phase = 0,
    data_hist = h1_data_combined,
    fitted_hists = fit_hists,
    outputpath = fit_pull_outputpath,
    component_labels = proc_labels,
    title = "Combined Informed Fit all phases",
    xlabel = "Sum Energy (keV)",
    ylabel = "Counts / 100 keV",
    residual_ylabel = "Data - Fit",
    pull_ylabel = "Pull",
    main_limits = (0, 4000, 0, nothing),
    pull_limits = (-4.5, 4.5),
    blinded_roi = (2700, 3000),
    figure_size = (2800, 1900),
    fontsize = 40,
    colors = plot_style_colors,
    total_fit_color = :black,
    total_fit_linewidth = 7,
    component_linewidth = 4,
    data_markersize = 14,
    data_whiskerwidth = 9,
    residual_band_color = (:gray65, 0.35),
    pull_band_color = (:gray75, 0.4),
    legend_position = :rt,
    show_chi2 = true,
    chi2_nparams = length(μ_vec),
    chi2_text_pos = (0.70, 0.88),
    show_pull_stats = true,
)



# plot_fit(
#     phase = 0,
#     data_hist = h1_data_combined,
#     fitted_hists = fit_hists,
#     fit_params = post_mean,
#     component_count_uncertainties = μ_vec,
#     show_component_uncertainties = true,
#     title = "Combined Informed Fit all phases",
#     outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_informed_combined.png"),
#     component_labels = proc_labels,
#     main_limits = (0, 180, 0, nothing),
#     ratio_limits = (0, 180, 0., 2.0),
#     # blinded_roi = (2700, 3000),
#     chi2_nparams = length(μ_vec),
#     include_component_counts = true,
#     figure_size = (2000, 1200),
#     chi2_text_pos = (0.1, 0.7),
#     fontsize = 42,
#     legend_tellheight=true,
#     legend_width = 650
# )
    display(p2)
    display(p1)
end
# Print summary of results



