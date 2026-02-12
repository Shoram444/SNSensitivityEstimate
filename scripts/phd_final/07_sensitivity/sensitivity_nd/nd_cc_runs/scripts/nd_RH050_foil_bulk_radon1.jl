using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Metaheuristics

println("loaded pkgs")


## Helper functions
include(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/helper_functions.jl"))


analysisDict = Dict(
    :signal => "RH050_foil_bulk",
    :radon_tag => 1,
    :hours => 22,
    :side => %SIDE,
)
# analysisDict = Dict(
#     :signal => "Nnubb500keV_foil_bulk",
#     :radon_tag => 1,
#     :hours => 0.001,
#     :side => "sameSide",
# )

if analysisDict[:side] == "both"
    vars = [
        "phi", 
        "sumE", 
        "dy", 
        "dz",
        "lPint",
        "lPext",
        ]

    bin_bounds = (
        phi = (0,180),
        sumE = (300, 3000),
        dy = (0, 150),
        dz = (0, 150),
        lPint = (0, 10),
        lPext = (0, 50),
    )
else 
    vars = [
        "phi", 
        "sumE", 
        "dy", 
        "dz",
        "sameSide",
        "lPint",
        "lPext",
        ]

    bin_bounds = (
        phi = (0,180),
        sumE = (300, 3000),
        dy = (0, 150),
        dz = (0, 150),
        sameSide = (0, 1),
        lPint = (0, 10),
        lPext = (0, 50),
    )
end

processes = load_ndim_processes(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), bin_bounds, vars)

data_dir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/neutrons_jan_2026/")

include(joinpath(data_dir, "read_neutrons_ND.jl"))
neutron_processes = load_neutron_processes_ND(data_dir, vars, bin_bounds)
for p in neutron_processes
    set_activity!(p, p.activity / 3)
    println("Neutron process: ", p.isotopeName, " nTotalSim: ", p.nTotalSim, " activity: ", p.activity)
end


# signal_name = "RH050_foil_bulk"
signal_name = analysisDict[:signal]
signal = get_process(signal_name, processes) |> first
# signal = get_process("RH037_foil_bulk", processes) |> first
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
bkgs = [    
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


background = [get_process(bkg, processes) |> first for bkg in bkgs]
background = vcat(background, neutron_processes)

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)
# set_activity!(background[7], 0.15/1e3) # radon to 0.15mBq/kg
# set_activity!(background[7], 2/1e3) # radon to 2mBq/kg
radon_tag = analysisDict[:radon_tag]
if radon_tag == 1
    set_activity!(background[7], 150/1e6) # radon to 150 uBq/kg
elseif radon_tag == 2
    set_activity!(background[7], 2/1e3) # radon to 2 mBq/kg
elseif radon_tag == 3
    set_activity!(background[7], 0.6/1e3) # radon to 0.6 mBq/kg
end

const α_ = 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


function filter_process!(process::DataProcessND, var_name::Symbol, range::Tuple{Real, Real})
    process.data = filter(x -> getproperty(x, var_name) >= range[1] && getproperty(x, var_name) <= range[2], process.data)
    @info "filtered process $(process.isotopeName) on $var_name in range $range, new nEvents: $(length(process.data))"
    return process;
end

if (analysisDict[:side] == "sameSide") 
    side_range = (1., 1.) 
    filter_process!(signal, :sameSide, side_range);
    map( x-> filter_process!(x, :sameSide, side_range), background);

elseif (analysisDict[:side] == "oppositeSide") 
    side_range = (0., 0.) 
    filter_process!(signal, :sameSide, side_range);
    map( x-> filter_process!(x, :sameSide, side_range), background);
end



searchRange = make_stepRange(signal)


# manual bounds
if analysisDict[:side] == "both"
    lb = float.([0, 140, 300, 2000, 0, 50, 0, 50, 0, 0, 0, 50 ])
    ub = float.([80, 180, 1600, 2800, 10, 150, 10, 150, 0, 10, 10, 50])
elseif analysisDict[:side] == "sameSide"
    lb = float.([0, 0, 300, 2000, 0, 50, 0, 50, 0,2, 0, 0, 0, 50 ])
    ub = float.([180, 180, 1600, 2800, 10, 150, 10, 150,  0,2,0, 10, 10, 50])
else 
    lb = float.([0, 0, 300, 2000, 0, 50, 0, 50, 0,1, 0, 0, 0, 50 ])
    ub = float.([180, 180, 1600, 2800, 10, 150, 10, 150,  0,1,0, 10, 10, 50])    
end

PROCESSES = vcat(signal, background)


function prob(x::Vector{Float64})
    # xs = normalized_to_roi_vector(collect(x))
    return - get_s_to_b(SNparams, α_, PROCESSES, x; approximate="formula")
end



options = Options(;
    x_tol = 1e-2,
    f_tol = 1e-4,
    f_tol_rel = 1e-5,
    f_tol_abs = 1e-5,
    time_limit = analysisDict[:hours]*60*60.0, #60*60*12.0,
    parallel_evaluation = false,
    verbose = true,
)


box_bounds = boxconstraints(lb = lb, ub = ub)

algo = ECA(N=20, K=3, ;options)

if analysisDict[:side] == "both"
    x0 = ( float.( [0, 180, 400, 2400, 0, 120, 0, 120, 0, 4, 1, 50] ) )
    set_user_solutions!(algo, x0, prob)
elseif analysisDict[:side] == "sameSide"
    x0 = ( float.( [0, 180, 400, 2400, 0, 120, 0, 120, 0,2, 0,3, 1, 50] ) )
    set_user_solutions!(algo, x0, prob)
elseif analysisDict[:side] == "oppositeSide"
    x0 = ( float.( [0, 180, 400, 2400, 0, 120, 0, 120, 0,1, 0,3, 1, 50] ) )
    set_user_solutions!(algo, x0, prob)
end


result = Metaheuristics.optimize(prob, box_bounds, algo)
@show minimum(result)
@show res=  minimizer(result)

best = get_best_ROI_ND(res, signal)
best_sens = get_sensitivityND(SNparams, α_, PROCESSES, best; approximate="table")
print(best_sens)


