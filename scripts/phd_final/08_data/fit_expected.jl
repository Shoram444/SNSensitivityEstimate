using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using UnROOT

println("loaded pkgs")

include(datadir("sims/final_phd/data_sims_mul95_alpha_veto_curved_pcd/load_nd_sim_files.jl"))
include(scriptsdir("phd_final/08_data/helper_functions.jl"))


data_dir = datadir("data/final_phd/data_no_gamma_alpha_veto_nonan_fixed_us_window")
f1 = ROOTFile(joinpath(data_dir, "phase1_data.root"))
d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

f2 = ROOTFile(joinpath(data_dir, "phase2_data.root"))
d2 = LazyTree(f2, "tree", keys(f2["tree"])) |> DataFrame

f3 = ROOTFile(joinpath(data_dir, "phase3_data.root"))
d3 = LazyTree(f3, "tree", keys(f3["tree"])) |> DataFrame

fit_results_exp = []

for p =1:3
    phase = p
    var_data = :sumE
    var_simu = :sumEsimu
    binning = (0:100:4000)
    fwhm = 0.12
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

    e1_min, e1_max = data_cuts[:e1]

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
    simu_p3_processes_roi_1 = load_sim_processes_ND(p1_dir, simu_varNames_1, simu_roi_1; timeMeas = p3_duration_seconds)


    d1_roi_1 = filter_data!(d1, data_roi_1)
    d2_roi_1 = filter_data!(d2, data_roi_1)
    d3_roi_1 = filter_data!(d3, data_roi_1)
    

    h1_data_p1 = Hist1D(getproperty(d1_roi_1, var_data); binedges = binning) #* inv(p1_duration_days)
    h1_data_p2 = Hist1D(getproperty(d2_roi_1, var_data); binedges = binning) #* inv(p2_duration_days)
    h1_data_p3 = Hist1D(getproperty(d3_roi_1, var_data); binedges = binning) #* inv(p3_duration_days)

    # h_data_tot = h1_data_p1 + h1_data_p2 + h1_data_p3
    # h_data_tot = (h1_data_p1 + h1_data_p2) * 0.5
    h_data_tot = h1_data_p2


    p1_K40_process = get_process("K40_foil_bulk", simu_p1_processes_roi_1) |> first
    p1_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p1_processes_roi_1) |> first
    p1_bb_process = get_process("bb_foil_bulk", simu_p1_processes_roi_1) |> first
    p1_radon_process = get_process("Bi214_wire_surface", simu_p1_processes_roi_1) |> first
    # p1_K40_PMT_process = get_process("K40_PMT_glass_bulk", simu_p1_processes_roi_1) |> first

    p2_K40_process = get_process("K40_foil_bulk", simu_p2_processes_roi_1) |> first
    p2_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p2_processes_roi_1) |> first
    p2_bb_process = get_process("bb_foil_bulk", simu_p2_processes_roi_1) |> first
    p2_radon_process = get_process("Bi214_wire_surface", simu_p2_processes_roi_1) |> first


    p3_K40_process = get_process("K40_foil_bulk", simu_p3_processes_roi_1) |> first
    p3_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p3_processes_roi_1) |> first
    p3_bb_process = get_process("bb_foil_bulk", simu_p3_processes_roi_1) |> first
    p3_radon_process = get_process("Bi214_wire_surface", simu_p3_processes_roi_1) |> first

    set_activity!(p1_radon_process, 0.1011)
    set_activity!(p2_radon_process, 0.0677)
    set_activity!(p3_radon_process, 3.36e-3)
    set_activity!(p1_Pa234m_process, 2.85e-3) # 2.85 vs 17.3
    set_activity!(p2_Pa234m_process, 2.85e-3) # 2.85 vs 17.3
    set_activity!(p3_Pa234m_process, 2.85e-3) # 2.85 vs 17.3
    

    simu_h_p1_K40 = get_expected_simu_counts(p1_K40_process, var_simu, fwhm, bins = binning, e_range = (e1_min, e1_max)) *1.0
    simu_h_p1_Pa234m = get_expected_simu_counts(p1_Pa234m_process, var_simu, fwhm, bins = binning, e_range = (e1_min, e1_max)) *1.0
    simu_h_p1_bb = get_expected_simu_counts(p1_bb_process, var_simu, fwhm, bins = binning, e_range = (e1_min, e1_max)) *1.0
    # simu_h_p1_K40_PMT = get_expected_simu_counts(p1_K40_PMT_process, var_simu, fwhm, bins = binning, e_range = (e1_min, e1_max)) * 1.0
    simu_h_p1_radon = get_expected_simu_counts(p1_radon_process, var_simu, fwhm, bins = binning, e_range = (e1_min, e1_max)) *1.0

    simu_h_p2_K40 = get_expected_simu_counts(p2_K40_process, var_simu, fwhm, bins = binning, timeMeas = p2_duration_seconds, e_range = (e1_min, e1_max)) *1.0
    simu_h_p2_Pa234m = get_expected_simu_counts(p2_Pa234m_process, var_simu, fwhm, bins = binning, timeMeas = p2_duration_seconds, e_range = (e1_min, e1_max)) *1.0
    simu_h_p2_bb = get_expected_simu_counts(p2_bb_process, var_simu, fwhm, bins = binning, timeMeas = p2_duration_seconds, e_range = (e1_min, e1_max)) *1.0
    # simu_h_p2_K40_PMT = get_expected_simu_counts(p2_K40_PMT_process, var_simu, fwhm, bins = binning, timeMeas = p2_duration_seconds, e_range = (e1_min, e1_max)) * 1.0
    simu_h_p2_radon = get_expected_simu_counts(p2_radon_process, var_simu, fwhm, bins = binning, timeMeas = p2_duration_seconds, e_range = (e1_min, e1_max)) *1.0

    simu_h_p3_K40 = get_expected_simu_counts(p3_K40_process, var_simu, fwhm, bins = binning, timeMeas = p3_duration_seconds, e_range = (e1_min, e1_max)) *1.0
    simu_h_p3_Pa234m = get_expected_simu_counts(p3_Pa234m_process, var_simu, fwhm, bins = binning, timeMeas = p3_duration_seconds, e_range = (e1_min, e1_max)) *1.0
    simu_h_p3_bb = get_expected_simu_counts(p3_bb_process, var_simu, fwhm, bins = binning, timeMeas = p3_duration_seconds, e_range = (e1_min, e1_max)) *1.0
    # simu_h_p3_K40_PMT = get_expected_simu_counts(p3_K40_PMT_process, var_simu, fwhm, bins = binning, timeMeas = p3_duration_seconds, e_range = (e1_min, e1_max)) * 1.0
    simu_h_p3_radon = get_expected_simu_counts(p3_radon_process, var_simu, fwhm, bins = binning, timeMeas = p3_duration_seconds, e_range = (e1_min, e1_max)) *1.0

    radon_level1 = p1_radon_process.activity * 1000 # in mBq/m3
    radon_level2 = p2_radon_process.activity * 1000 # in mBq/m3
    radon_level3 = p3_radon_process.activity * 1000 # in mBq/m3
    Pa_level = p1_Pa234m_process.activity * 1000 # in mBq/kg

    # h_flat = Hist1D(rand(Uniform(500, 3000),1000); binedges = binning)

    # tot_duration = p1_duration_seconds + p2_duration_seconds + p3_duration_seconds


    data_hist = if phase == 1
        h1_data_p1
    elseif phase == 2
        h1_data_p2
    elseif phase == 3
        h1_data_p3
    end

    bb = if phase == 1
        simu_h_p1_bb
    elseif phase == 2
        simu_h_p2_bb
    elseif phase == 3
        simu_h_p3_bb
    end

    K40 = if phase == 1
        simu_h_p1_K40
    elseif phase == 2
        simu_h_p2_K40
    elseif phase == 3
        simu_h_p3_K40
    end

    Pa = if phase == 1
        simu_h_p1_Pa234m
    elseif phase == 2
        simu_h_p2_Pa234m
    elseif phase == 3
        simu_h_p3_Pa234m
    end

    radon = if phase == 1
        simu_h_p1_radon
    elseif phase == 2
        simu_h_p2_radon
    elseif phase == 3
        simu_h_p3_radon
    end

    eff_bb = if phase == 1
        length(p1_bb_process.data) / p1_bb_process.nTotalSim
    elseif phase == 2
        length(p2_bb_process.data) / p2_bb_process.nTotalSim
    elseif phase == 3
        length(p3_bb_process.data) / p3_bb_process.nTotalSim
    end

    simu_h_total = K40 + Pa + bb + radon #+ h_flat #+ K40_PMT
    simu_combined = [bb, K40, Pa, radon] #+ h_flat #+ K40_PMT

    labels = [
        rich("2νββ"),
        rich(superscript("40"), "K"),
        rich(superscript("234m"), "Pa"),
        rich("radon"),
    ]

    activities = (
        bb = p1_bb_process.activity,
        K40 = p1_K40_process.activity,
        Pa = Pa_level / 1000,
        radon = radon_level1 / 1000,
    )

    res = fit_result(phase, data_hist, simu_combined, activities, labels)
    push!(fit_results_exp, res)

    plot_fit(
        phase = phase,
        data_hist = data_hist,
        fitted_hists = simu_combined,
        fit_params = nothing,
        outputpath = scriptsdir("phd_final/08_data/figs/expected_phase_$(phase).png"),
        component_labels = labels,
        main_limits = (0, 4000, 0, nothing),
        ratio_limits = (0, 4000, 0.2, 1.8),
        blinded_roi = (2700, 3000),
        chi2_nparams = length(simu_combined),
        include_component_counts = true,
    )

end

data_combined = sum([p.data_hist for p in fit_results_exp])

fit_p1 = fit_results_exp[1].fitted_hists
fit_p2 = fit_results_exp[2].fitted_hists
fit_p3 = fit_results_exp[3].fitted_hists

fit_combined = sum([fit_p1, fit_p2, fit_p3])
fit_total = sum([sum(fit_p1), sum(fit_p2), sum(fit_p3)])

f_full = plot_fit(
    phase = 0,
    data_hist = data_combined,
    fitted_hists = fit_combined,
    fit_params = nothing,
    outputpath = scriptsdir("phd_final/08_data/figs/bayes_fit_combined_expected.png"),
    component_labels = proc_labels,
    main_limits = (0, 4000, 0, nothing),
    ratio_limits = (0, 4000, 0.2, 1.8),
    blinded_roi = (2700, 3000),
    chi2_nparams = length(fit_combined),
    title = "Combined Expected phase 1+2+3",
    include_component_counts = true,
)