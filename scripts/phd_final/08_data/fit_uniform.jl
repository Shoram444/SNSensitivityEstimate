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

begin
data_dir = datadir("data/final_phd/data_no_gamma_alpha_veto_nonan_fixed_us_window")
f1 = ROOTFile(joinpath(data_dir, "phase1_data.root"))
d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

f2 = ROOTFile(joinpath(data_dir, "phase2_data.root"))
d2 = LazyTree(f2, "tree", keys(f2["tree"])) |> DataFrame

f3 = ROOTFile(joinpath(data_dir, "phase3_data.root"))
d3 = LazyTree(f3, "tree", keys(f3["tree"])) |> DataFrame

var_data = :sumE
var_simu = :sumEsimu
binning = (0:50:3000)
fwhm = 0.15
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
end

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
    std_post = std(samples)

    mm = map(measurement, post_mean, std_post)

    μ_K40 = mm.μ_K40
    μ_Pa = mm.μ_Pa
    μ_bb = mm.μ_bb
    μ_radon = mm.μ_radon
    μ_internal = mm.μ_internal
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
            push!(hists, Measurements.value(μ[i]) * bkg_hists_normed[i])
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
            μ_vec,
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
        component_count_uncertainties = μ_vec,
        show_component_uncertainties = true,
        outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_phase$(phase)_free.png"),
        component_labels = proc_labels,
        main_limits = (0, 4000, 0, nothing),
        ratio_limits = (0, 4000, 0., 2.0),
        blinded_roi = (2700, 3000),
        chi2_nparams = length(μ_vec),
        include_component_counts = true,
        figure_size = (2400, 1250),
        chi2_text_pos = (0.77, 0.7),
        fontsize = 42,
        legend_tellheight=true,
        legend_width = 650
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

n_tot = fit_results_free[1].fit_n + fit_results_free[2].fit_n + fit_results_free[3].fit_n
combined_component_stds_free = Measurements.uncertainty.(n_tot)


f_full = plot_fit(
    phase = 0,
    data_hist = full_data,
    fitted_hists = fit_combined,
    fit_params = nothing,
    component_count_uncertainties = combined_component_stds_free,
    show_component_uncertainties = true,
    outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_combined_free_lines.png"),
    component_labels = proc_labels,
    main_limits = (0, 4000, 0, nothing),
    ratio_limits = (0, 4000, 0., 2.),
    blinded_roi = (2700, 3000),
    chi2_nparams = length(fit_combined),
    title = "Combined unconstrained fit phase 1+2+3",
    include_component_counts = true,
    stacked = false,
    filled = false,
    figure_size = (2400, 1250),
    chi2_text_pos = (0.77, 0.7),
    fontsize = 42,
    legend_tellheight=true,
    legend_width = 650
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
    margin = 5Plots.mm,
    figure_title = "Posterior samples of μ_bb vs μ_radon for phase 1",
    thickness_scaling = 1.3,
    size = (1200, 800),
    tickfontsize = 16,
    labelfontsize = 18,
    xrotation = 45,
)
savefig(p1, scriptsdir("phd_final/08_data/figs/unconstrained_parameter_cor_phase1_bb_radon.png"))

p2 = StatsPlots.plot(
    samples_by_phase[1],
    vsel = (:μ_Pa, :μ_K40),
    margin = 5Plots.mm,
    figure_title = "Posterior samples of μ_Pa vs μ_K40 for phase 1",
    thickness_scaling = 1.3,
    size = (1200, 800),
    tickfontsize = 16,
    labelfontsize = 18,
    xrotation = 45,
)
savefig(p2, scriptsdir("phd_final/08_data/figs/unconstrained_parameter_cor_phase1_Pa_K40.png"))



### activities
function get_activity_from_exp_n(a_exp, n_exp, n_fit)
    a_exp * n_fit / n_exp
end

total_duration = p1_duration_seconds + p2_duration_seconds + p3_duration_seconds
radon_a_exp = (p1_radon_process.activity * p1_duration_seconds + p2_radon_process.activity * p2_duration_seconds + p3_radon_process.activity * p3_duration_seconds) / total_duration
n_exp_radon = fit_results_exp[1].exp_n[4] + fit_results_exp[2].exp_n[4] + fit_results_exp[3].exp_n[4]
n_fit_radon = fit_results_free[1].fit_n[5] + fit_results_free[2].fit_n[5] + fit_results_free[3].fit_n[5]
a_radon_combined = get_activity_from_exp_n(radon_a_exp, n_exp_radon, n_fit_radon)

bb_a_exp = p1_bb_process.activity
n_exp_bb = fit_results_exp[1].exp_n[1] + fit_results_exp[2].exp_n[1] + fit_results_exp[3].exp_n[1]
n_fit_bb = fit_results_free[1].fit_n[1] + fit_results_free[2].fit_n[1] + fit_results_free[3].fit_n[1]
a_bb_combined = get_activity_from_exp_n(bb_a_exp, n_exp_bb, n_fit_bb)

K40_a_exp = p1_K40_process.activity
n_exp_K40 = fit_results_exp[1].exp_n[2] + fit_results_exp[2].exp_n[2] + fit_results_exp[3].exp_n[2]
n_fit_K40 = fit_results_free[1].fit_n[2] + fit_results_free[2].fit_n[2] + fit_results_free[3].fit_n[2]
a_K40_combined = get_activity_from_exp_n(K40_a_exp, n_exp_K40, n_fit_K40)

Pa_a_exp = 2.85e-3
n_exp_Pa = fit_results_exp[1].exp_n[3] + fit_results_exp[2].exp_n[3] + fit_results_exp[3].exp_n[3]
n_fit_Pa = fit_results_free[1].fit_n[3] + fit_results_free[2].fit_n[3] + fit_results_free[3].fit_n[3]
a_Pa_combined = get_activity_from_exp_n(Pa_a_exp, n_exp_Pa, n_fit_Pa)



using PrettyTables
d = CSV.File(scriptsdir("phd_final/08_data/data_mc_numbers.csv")) |> DataFrame

@transform! d :a_unconstrained = measurement.(d.a_unconstrained)
@transform! d :a_constrained = measurement.(d.a_constrained)

pretty_table(d, backend = Val(:latex), header = ["Process", "Expected activity (Bq)", "Unconstrained fit activity (Bq)", "Constrained fit activity (Bq)"], alignment = :l, title = "Comparison of expected and fitted activities for different processes. The unconstrained fit allows all components to vary freely, while the constrained fit incorporates prior information on the activities. The values are given with their uncertainties.", title_alignment = :c)

@transform! d :diff_unconstrained = :a_unconstrained .- :a_exp
@transform! d :diff_constrained = :a_constrained .- :a_exp
