using DrWatson
@quickactivate "SNSensitivityEstimate"

using SNSensitivityEstimate
using Random, LinearAlgebra, Statistics, Distributions, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
using DataFramesMeta, CSV

global t0 = time()

# Configuration for sterile neutrino mass scan
sterile_masses = 800:100:1500  # 800, 900, 1000, 1100, 1200, 1300, 1400, 1500
n_samples_per_mass = 300

# Standard analysis parameters
bin_width = 100
Bin_low, Bin_high = 300, 3500
mode = "sumE"
prior_value = 1
radon_tag = 1

# Include ROI definitions
include(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/results/20260123_results/best_rois.jl"))

# Background processes (defined once, reused for all masses)
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

# Helper function for background histograms
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

        if p.isotopeName in internals && integral(h) > 0
            h_internal = h_internal + h
        elseif p.isotopeName in bb && integral(h) > 0
            h_bb = h_bb + h
        elseif p.isotopeName in radon && integral(h) > 0
            h_radon = h_radon + h
        elseif p.isotopeName in Bi210 && integral(h) > 0
            h_Bi210 = h_Bi210 + h
        elseif p.isotopeName in detectors && integral(h) > 0
            h_detector = h_detector + h
        elseif integral(h) > 0
            h_external = h_external + h
        end
    end
    return vcat(h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external)
end

# Load processes once
@info "Loading data processes..."
all_processes = load_data_processes(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), mode, fwhm=0.0, roi = ())

neutron_data_dir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/neutrons_jan_2026/")
include(joinpath(neutron_data_dir, "read_neutrons_1D.jl"))
neutron_processes = load_neutron_process1D(neutron_data_dir, mode)
for p in neutron_processes
    set_activity!(p, p.activity / 3)
end

@info "Data processes loaded."

# Dictionary to store results for each mass
results_by_mass = Dict()

# Loop over each sterile mass
for mass in sterile_masses
    println("\n" * "="^80)
    println("Processing sterile neutrino mass: $(mass) keV")
    println("="^80)

    sname = "sterile$(lpad(mass, 4, '0'))_foil_bulk"
    println("Signal process name: $sname")
    
    signal_name = sname
    
    # Set up ROI based on signal type
    roi_sumE = (Bin_low, Bin_high)
    roi = Dict(:sumE => roi_sumE)
    
    # Get signal process
    try
        global signal = get_process(signal_name, all_processes) |> first
    catch e
        @warn "Failed to load signal $signal_name: $e"
        continue
    end
    
    # Get background processes
    global background = [get_process(b, all_processes) |> first for b in backgrounds]
    background = vcat(background, neutron_processes)
    
    # Set radon activity
    if radon_tag == 1
        set_activity!(background[7], 150/1e6)
    elseif radon_tag == 2
        set_activity!(background[7], 2/1e3)
    elseif radon_tag == 3
        set_activity!(background[7], 0.6/1e3)
    end
    
    set_signal!(background[1], false)
    
    # Set bins
    for b in background
        set_bins!(b, Bin_low:bin_width:Bin_high+bin_width)
    end
    set_bins!(signal, Bin_low:bin_width:Bin_high+bin_width)
    
    # Get background histograms
    bkg_hist = get_combined_bkg_hists(background, Bin_low, Bin_high, bin_width)
    bkg_hist_normed = normalize.(bkg_hist, width = true)
    signal_hist_normed = normalize(restrict(get_bkg_counts_1D(signal), Bin_low, Bin_high+bin_width), width = true)
    
    # Define prior
    prior = NamedTupleDist(
        As = Uniform(1e-20, prior_value),
        Ab = [Uniform(1e-20, 1) for _ in 1:length(bkg_hist)] 
    )
    
    # Run MCMC sampling
    t_halfs = Float64[]
    println("Running $(n_samples_per_mass) MCMC samples...")
    
    for i in 1:n_samples_per_mass
        t1 = time()
        try 
            sens = get_sens_bayes_uniform(bkg_hist, signal, prior; 
                ROI_a = Bin_low, ROI_b = Bin_high, 
                nsteps = 10^4, nchains = 4, smallest_interval = true)
            push!(t_halfs, sens)
            
            if mod(i, 5) == 0
                @info "Sample $i/$n_samples_per_mass completed ($(time() - t1)s)"
            end
        catch e
            @warn "Sample $i failed: $e"
            continue
        end
    end
    
    # Store results in a DataFrame for this mass
    df_mass = DataFrame(
        mass_keV = fill(mass, length(t_halfs)),
        sensitivity_yr = t_halfs
    )
    results_by_mass[mass] = df_mass
    
    println("Completed $(length(t_halfs)) samples for mass $(mass) keV")
    println("Mean sensitivity: $(mean(t_halfs)) years")
    println("Std sensitivity: $(std(t_halfs)) years")
    
    # Save individual results for this mass
    save_name = "bayes_sterile$(lpad(mass, 4, '0'))_foil_bulk_radon$(radon_tag)_$(n_samples_per_mass)samples"
    CSV.write(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/sensitivities_$(save_name)_$(rand(1:1000000)).csv"), df_mass)
end

# Combine all results into a single DataFrame
println("\n" * "="^80)
println("Combining results...")
println("="^80)

df_combined = vcat(collect(values(results_by_mass))...)

# Add summary statistics
println("\nSummary Statistics:")
println(combine(groupby(df_combined, :mass_keV), 
    :sensitivity_yr => mean => :mean_sensitivity,
    :sensitivity_yr => std => :std_sensitivity,
    :sensitivity_yr => length => :n_samples
))

# Save combined results
save_name_combined = "bayes_sterile_scan_radon$(radon_tag)_$(n_samples_per_mass)samples_per_mass"
CSV.write(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/sensitivities_$(save_name_combined)_$(rand(1:1000000)).csv"), df_combined)

println("\nResults saved.")


# read saved results into one df
d = []
for f in readdir(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/"))
    if occursin("foil_bulk_radon$(radon_tag)_$(n_samples_per_mass)samples", f)
        df = CSV.read(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/$(f)"), DataFrame)
        push!(d, df)
    end
end
df_all = vcat(d...)
df_all = @chain df_all begin
    @groupby :mass_keV
    @combine :mean_sensitivity = mean(:sensitivity_yr)
    @orderby :mass_keV
end

df_all.mean_sensitivity