using DrWatson
@quickactivate "SNSensitivityEstimate"


using SNSensitivityEstimate
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
using DataFramesMeta, CSV

global t0 = time()



analysisDict = Dict(
    :signal => "bb0nuM2_foil_bulk",
    :bining => (300, 3500),
    :bin_width => 100,
    :mode => "sumE",
    :prior => 1e-3,# 1e-4 0nu, 1e-4 RH, 1e-3 M1, 1e-2 M2
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
# background = vcat(backgrounds, neutron_processes)


# declare which process is signal
signal = get_process(analysisDict[:signal], all_processes) |> first

# declare background processes
background = [get_process(b, all_processes) |> first for b in backgrounds]
background = vcat(background, neutron_processes)

if analysisDict[:radon_tag] == 1
    set_activity!(background[7], 150/1e6) # radon to 150 uBq/kg
elseif analysisDict[:radon_tag] == 2
    set_activity!(background[7], 2/1e3) # radon to 2 mBq/kg
elseif analysisDict[:radon_tag] == 3
    set_activity!(background[7], 0.6/1e3) # radon to 0.6 mBq/kg
end


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
# f2(pars::NamedTuple{(:As, :Ab)}, x::Real) = f_uniform_bkg(pars, x, signal_hist_normed, bkg_hist_normed)

# uninformed prior for each activity
prior = NamedTupleDist(
    As = Uniform(1e-20, analysisDict[:prior]), # 1e-4 0nu, 1e-4 RH, 1e-3 M1, 1e-2 M2
    Ab = [Uniform(1e-20,1) for _ in 1:length(bkg_hist)] 
)   

t_halfs = Float64[]
while(time() - t0 < 3600*12) # do this for n hours
# for _ in 1:100 # do this for n hours
# while(time() - t0 < 3600*12) # do this for n hours
# for _ in 1:1 # do this for n hours
    # GC.gc()
    t1 = time()
    try 
        sens = get_sens_bayes_uniform(bkg_hist, signal, prior; ROI_a = Bin_low, ROI_b = Bin_high, nsteps = 5*10^4, nchains = 4)
        println("time to fit, t = $(time() - t1) s")
        println(sens)
        push!(t_halfs, sens)
        # CSV.write(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/nu0bb/sensitivities_$(save_name)_$(rand(1:1000000)).csv"), DataFrame(thalf= t_halfs))
    catch
        @warn "failed fit" 
        continue
    end
end

# histogram(t_halfs, bins = 20, xlabel = "Sensitivity (years)", ylabel = "Frequency", title = "Bayesian Sensitivity Estimates for $(analysisDict[:signal]) with prior $(analysisDict[:prior]) and ROI [$(Bin_low), $(Bin_high)]")
save_name = savename(analysisDict)
CSV.write(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/nu0bb/sensitivities_$(save_name)_$(rand(1:1000000)).csv"), DataFrame(thalf= t_halfs))
