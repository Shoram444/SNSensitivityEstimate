using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using UnROOT
using BAT, SpecialFunctions, DensityInterface, IntervalSets, ValueShapes
using StatsPlots

println("loaded pkgs")

include(datadir("sims/final_phd/data_sims_mul95_alpha_veto_curved_pcd/load_nd_sim_files.jl"))
include(scriptsdir("phd_final/08_data/bayes_fit_models.jl"))

function filter_data!(d, roi)
    d_new = deepcopy(d)
    for key in keys(roi)
        lo, hi = roi[key]
        filter!(row -> lo <= row[key] <= hi, d_new)
    end
    return d_new
end

function smear_energy(E::Real, fwhm::Real) 
    if(fwhm == 0)
        return E
    end
    sigma = fwhm/2.355 * E * sqrt(1000/E)
    return rand(Normal(E, sigma))
end

function get_smeared_energy(p, var, fwhm)
    if var == :simulatedEnergy1 || var == :simulatedEnergy2
        return smear_energy.(getproperty.(p.data, var), fwhm)
    elseif var == :sumEsimu
        return smear_energy.(getproperty.(p.data, :simulatedEnergy1), fwhm) + smear_energy.(getproperty.(p.data, :simulatedEnergy2), fwhm)
    else
        return getproperty.(p.data, var)
    end
end

function get_expected_simu_counts(p, var, fwhm; bins = (0:100:4000), timeMeas = p.timeMeas)
    smeared_energy = fwhm == 0 ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm)
    h = normalize(Hist1D(smeared_energy; binedges = bins), width = false)
    n_exp = p.activity * timeMeas * p.amount * (length(p.data) / p.nTotalSim)
    return h * n_exp
end

function get_simu_histo(p, var, fwhm; bins = (0:100:4000))
    smeared_energy = fwhm == 0 ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm)
    h = Hist1D(smeared_energy; binedges = bins)
    return h 
end

function get_chi2_ndf(data, fit, nparams)
    bincount_data = bincounts(data)
    bincount_fit = bincounts(fit)

    binerrors_data = binerrors(data)
    binerrors_fit = binerrors(fit)

    chi2 = 0.0
    ndf = 0

    for i in eachindex(bincount_data)
        O = bincount_data[i]
        E = bincount_fit[i]

        s = binerrors_data[i]^2 + binerrors_fit[i]^2

        if E > 0          # avoid division by zero
            chi2 += (O - E)^2 / s
            ndf += 1
        end
    end

    ndf = ndf - nparams
    return chi2 / ndf, chi2, ndf
end

data_dir = datadir("data/final_phd/data_no_gamma_alpha_veto_no_nan")
f1 = ROOTFile(joinpath(data_dir, "phase1_data.root"))
d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

f2 = ROOTFile(joinpath(data_dir, "phase2_data.root"))
d2 = LazyTree(f2, "tree", keys(f2["tree"])) |> DataFrame

f3 = ROOTFile(joinpath(data_dir, "phase3_data.root"))
d3 = LazyTree(f3, "tree", keys(f3["tree"])) |> DataFrame

var_data = :sumE
var_simu = :sumEsimu
binning = (500:100:3000)
fwhm = 0.12
data_cuts = Dict(
    :sumE => (500, 3000),
    :e1 => (250, 5000),
    :e2 => (250, 5000),
    :dy => (0, 60),
    :dz => (0, 80),
    :trackLength1 => (300, 1100),
    :trackLength2 => (300, 1100),
    :deltaCaloTime => (0, 2.5),
    :phi => (0,180)
)



simu_roi_1 = (
    sumEsimu = data_cuts[:sumE],
    simulatedEnergy1 = data_cuts[:e1],
    simulatedEnergy2 = data_cuts[:e2],
    # dy = data_cuts[:dy],
    # dz = data_cuts[:dz],
    # trackLength1 = data_cuts[:trackLength1],
    # trackLength2 = data_cuts[:trackLength2],
    # deltaCaloTime = data_cuts[:deltaCaloTime],
    # phi = data_cuts[:phi]
)
simu_roi_2 = (
    sumE = data_cuts[:sumE],
    reconstructedEnergy1 = data_cuts[:e1],
    reconstructedEnergy2 = data_cuts[:e2],
    # dy = data_cuts[:dy],
    # dz = data_cuts[:dz],
    # trackLength1 = data_cuts[:trackLength1],
    # trackLength2 = data_cuts[:trackLength2],
    # deltaCaloTime = data_cuts[:deltaCaloTime],
    # phi = data_cuts[:phi]
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

p3_K40_process = get_process("K40_foil_bulk", simu_p2_processes_roi_1) |> first
p3_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p2_processes_roi_1) |> first
p3_bb_process = get_process("bb_foil_bulk", simu_p2_processes_roi_1) |> first
p3_radon_process = get_process("Bi214_wire_surface", simu_p1_processes_roi_1) |> first


simu_h_p1_K40 = get_simu_histo(p1_K40_process, var_simu, fwhm; bins = binning) 
simu_h_p1_Pa234m = get_simu_histo(p1_Pa234m_process, var_simu, fwhm; bins = binning)
simu_h_p1_bb = get_simu_histo(p1_bb_process, var_simu, fwhm; bins = binning)
simu_h_p1_radon = get_simu_histo(p1_radon_process, var_simu, fwhm; bins = binning)

simu_h_p2_K40 = get_simu_histo(p2_K40_process, var_simu, fwhm; bins = binning) 
simu_h_p2_Pa234m = get_simu_histo(p2_Pa234m_process, var_simu, fwhm; bins = binning)
simu_h_p2_bb = get_simu_histo(p2_bb_process, var_simu, fwhm; bins = binning)
simu_h_p2_radon = get_simu_histo(p2_radon_process, var_simu, fwhm; bins = binning)

simu_h_p3_K40 = get_simu_histo(p3_K40_process, var_simu, fwhm; bins = binning) 
simu_h_p3_Pa234m = get_simu_histo(p3_Pa234m_process, var_simu, fwhm; bins = binning)
simu_h_p3_bb = get_simu_histo(p3_bb_process, var_simu, fwhm; bins = binning)
simu_h_p3_radon = get_simu_histo(p3_radon_process, var_simu, fwhm; bins = binning)


h_flat_p2 = Hist1D(rand(Uniform(500, 3000),100000); binedges = binning)


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

bkg_hist = [
    simu_h_p2_K40,
    simu_h_p2_Pa234m,
    simu_h_p2_bb,
    simu_h_p1_radon,
    h_internal,
    h_Bi210,
    h_detector,
    h_external,
    # h_flat_p2,

]

data_hist = h1_data_p2
bkg_hist_normed = normalize.(bkg_hist, width = false)


##### Constraining some parameters

global C_K40 = p2_duration_seconds * p1_K40_process.amount
global C_Pa = p2_duration_seconds * p1_Pa234m_process.amount
global C_bb = p2_duration_seconds * p1_bb_process.amount
global C_radon = p2_duration_seconds * p1_radon_process.amount


global a_bb_fixed    = p1_bb_process.activity
global a_radon_fixed = 67.7 / 1000

function lognormal_from_mean_relerr(mean, relerr)
    σ_log = sqrt(log(1 + relerr^2))
    μ_log = log(mean) - 0.5 * σ_log^2
    return LogNormal(μ_log, σ_log)
end

function make_hist_likelihood_constrained(data_hist, bkg_hists)

    data_counts = bincounts(data_hist)
    bkg_counts = [bincounts(h) for h in bkg_hists]

    nbins = length(data_counts)

    logfuncdensity(function(p)

        # --- Compute μ from physics ---

        μ_K40   = C_K40   * p.a_K40   * p.eff_K40
        μ_Pa    = C_Pa    * p.a_Pa    * p.eff_Pa
        μ_bb    = C_bb    * a_bb_fixed    * p.eff_bb
        μ_radon = C_radon * a_radon_fixed * p.eff_radon
        # μ_flat  = p.μ_flat

        μ_internal = p.μ_internal
        μ_Bi210 = p.μ_Bi210
        μ_detector = p.μ_detector
        μ_external = p.μ_external

        μ_vec = (μ_K40, μ_Pa, μ_bb, μ_radon, μ_internal, μ_Bi210, μ_detector, μ_external)

        ll = 0.0

        for k in 1:nbins

            λk = 0.0

            for i in 1:8
                λk += μ_vec[i] * bkg_counts[i][k]
            end

            λk = max(λk, eps())
            ll += log_pdf_poisson(λk, data_counts[k])
        end

        return ll
    end)
end

using Distributions

# --- efficiencies (multiplicative uncertainties) ---
eff_K40_prior   = lognormal_from_mean_relerr(1e-4, 0.3)   # 30% relative uncertainty
eff_Pa_prior    = lognormal_from_mean_relerr(1e-3, 0.3)   
eff_radon_prior = lognormal_from_mean_relerr(1e-4, 0.3)  
eff_bb_prior = lognormal_from_mean_relerr(1e-1, 0.3)   # 30% relative uncertainty

# --- activities ---
a_K40_prior = Truncated(Normal(0.0587, 0.015), 0.0, Inf)
a_Pa_prior  = Truncated(Normal(0.01, 0.004), 0.0, Inf)
a_radon_prior = Truncated(Normal(67.7/1000, 0.001), 0.0, Inf)

# --- flat component ---
# μ_flat_prior = lognormal_from_mean_relerr(0.0, 0.4)
# μ_flat_prior = Uniform(0.0, 2*integral(data_hist))

μ_internal_prior = Uniform(0.0, 2*integral(data_hist))
μ_Bi210_prior = Uniform(0.0, 2*integral(data_hist))
μ_detector_prior = Uniform(0.0, 2*integral(data_hist))
μ_external_prior = Uniform(0.0, 2*integral(data_hist))


prior_constrained = NamedTupleDist(

    # --- K40 ---
    a_K40 = a_K40_prior,
    eff_K40 = eff_K40_prior,

    # --- Pa234m ---
    a_Pa = a_Pa_prior,
    eff_Pa = eff_Pa_prior,

    # --- bb ---
    eff_bb = eff_bb_prior,

    # --- radon ---
    a_radon = a_radon_prior,
    eff_radon = eff_radon_prior,

    # --- flat ---
    # μ_flat = μ_flat_prior

    # --- unconstrained components ---
    μ_internal = μ_internal_prior,
    μ_Bi210 = μ_Bi210_prior,
    μ_detector = μ_detector_prior,
    μ_external = μ_external_prior
)

likelihood_constrained = make_hist_likelihood_constrained(data_hist, bkg_hist_normed)
posterior_constrained = PosteriorMeasure(likelihood_constrained, prior_constrained)

burnin = MCMCMultiCycleBurnin(max_ncycles=100, nsteps_final=3000)
mcmcalgo = MetropolisHastings(
    weighting=RepetitionWeighting(),
    tuning=AdaptiveMHTuning()
)

samples, _ = bat_sample(
    posterior_constrained,
    MCMCSampling(
        mcalg=mcmcalgo,
        burnin=burnin,
        nsteps=5*10^5,
        nchains=4
    )
) 

binned_unshaped_samples, _ = bat_transform(Vector, samples)

Plots.plot(samples, vsel = (1:7),size = (2800, 2200), thickness_scaling = 1.3, margin= 6Plots.mm, dpi = 300)
savefig(current(), scriptsdir("phd_final/08_data/figs/bayes_fit_constrained_corner.png"))

post_mean = mean(samples)

a_K40   = post_mean.a_K40
eff_K40 = post_mean.eff_K40

a_Pa    = post_mean.a_Pa
eff_Pa  = post_mean.eff_Pa

eff_bb  = post_mean.eff_bb

a_radon   = a_radon_fixed
eff_radon = post_mean.eff_radon

# μ_flat = post_mean.μ_flat

μ_internal = post_mean.μ_internal
μ_Bi210 = post_mean.μ_Bi210
μ_detector = post_mean.μ_detector
μ_external = post_mean.μ_external

μ_vec = [
    C_K40   * a_K40   * eff_K40,
    C_Pa    * a_Pa    * eff_Pa,
    C_bb    * a_bb_fixed * eff_bb,
    C_radon * a_radon * eff_radon,
    # μ_flat  = μ_flat, --- IGNORE ---
    μ_internal,
    μ_Bi210,
    μ_detector,
    μ_external
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

activities = [
    a_K40,
    a_Pa,
    a_bb_fixed,
    a_radon,
    μ_internal,  # flat has no activity meaning
    μ_Bi210,
    μ_detector,
    μ_external,
    # NaN  # flat has no activity meaning
]

using CairoMakie

let
    colors = ColorSchemes.tab10 #[:orange, :green, :royalblue, :purple, :gray]
    proc_labels = [
        "K40",
        "Pa234m",
        "2nubb",
        "Bi214 (radon)",
        "internal",
        "Bi210",
        "detector",
        "external",
        # "flat"
    ]

    # --- figure layout ---
    f = Figure(size = (2000, 1400), fontsize = 28)

    ax = Axis(f[1,1],
        xlabel = "",
        ylabel = "Counts / 100 keV",
        title = "Phase 2",
        limits = (0, 4000, 0, nothing)
    )

    ax_ratio = Axis(f[2,1],
        xlabel = "Sum Energy [keV]",
        ylabel = "Data / Fit",
        limits = (0, 4000, 0.2, 1.8)
    )

    # -----------------------
    # DATA
    # -----------------------
    CairoMakie.scatter!(ax, data_hist,
        color = :black,
        markersize = 10,
        label = "Data"
    )

    CairoMakie.errorbars!(ax, data_hist,
        color = :black,
        whiskerwidth = 8
    )

    # -----------------------
    # COMPONENTS
    # -----------------------
    legend_entries = Any[]
    legend_labels = String[]

    activities = [a_K40, a_Pa, a_bb_fixed, a_radon, μ_internal, μ_Bi210, μ_detector, μ_external]
    efficiencies = [eff_K40, eff_Pa, eff_bb, eff_radon, NaN, NaN, NaN, NaN, NaN]

    for i in eachindex(fit_hists)

        CairoMakie.stephist!(
            ax,
            fit_hists[i],
            color = colors[i],
            linewidth = 3
        )

        n_fit = sum(bincounts(fit_hists[i]))

        if i ≤ 4
            lbl = string(
                proc_labels[i], "\n",
                "  n = ", round(n_fit, digits=1), "\n",
                "  a = ", round(activities[i], digits=5),
                (i == 3 ? " (fixed)" : ""),
                "\n",
                "  ε = ", round(efficiencies[i], digits=6),
                "\n"
            )
        else
            lbl = string(
                proc_labels[i], "\n",
                "  n = ", round(n_fit, digits=1),
                "\n"
            )
        end

        push!(legend_entries, LineElement(color = colors[i], linewidth=4))
        push!(legend_labels, lbl)
    end


    # -----------------------
    # TOTAL FIT
    # -----------------------
    CairoMakie.stephist!(
        ax,
        total_fit,
        color = :red,
        linewidth = 4,
        label = "Total Fit"
    )

    # -----------------------
    # BLINDED ROI
    # -----------------------
    CairoMakie.vspan!(ax, [2700], [3000], color = :black)
    CairoMakie.vspan!(ax_ratio, [2700], [3000], color = :black)

    # -----------------------
    # χ²
    # -----------------------
    chi_ndf = get_chi2_ndf(data_hist, total_fit, 7)

    chi_text = string(
        "χ² / ndf = ",
        round(chi_ndf[2], digits=2),
        " / ",
        chi_ndf[3],
        " = ",
        round(chi_ndf[1], digits=2)
    )

    text!(ax, 0.77, 0.9,
        text = chi_text,
        space = :relative,
        fontsize = 26
    )

    # -----------------------
    # RATIO PANEL
    # -----------------------
    x_min, x_max, x_step = binedges(data_hist)[1], binedges(data_hist)[end], binedges(data_hist)[2] - binedges(data_hist)[1]

    xs = collect(range(x_min + x_step/2, x_max - x_step/2, step = x_step))
    ys = bincounts(data_hist ./ total_fit)

    # centers = midpoints(bincounts(data_hist))

    CairoMakie.scatter!(ax_ratio, xs, ys, color=:black)
    CairoMakie.errorbars!(ax_ratio, xs, ys, binerrors(data_hist ./ total_fit), color=:black, whiskerwidth = 8)

    CairoMakie.hlines!(ax_ratio, [1.0], linestyle=:dash, color=:red)

    # -----------------------
    # LEGEND PANEL
    # -----------------------
    Legend(f[1,2], legend_entries, legend_labels,
        tellheight = false,
        patchsize = (30, 40, 100, 0),
        rowgap = 30,
    )

    rowsize!(f.layout, 2, Relative(0.25))

    save(
        scriptsdir("phd_final/08_data/figs/bayes_fit_phase2_pub.png"),
        f,
        px_per_unit = 3
    )

    f
end




################
