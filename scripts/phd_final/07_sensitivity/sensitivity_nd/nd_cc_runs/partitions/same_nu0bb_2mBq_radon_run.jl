using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Dates
now = Dates.now()
const signal_name = "bb0nu_foil_bulk"
const side = "sameSide" # "oppositeSide
const radon_tag = 2

println("loaded pkgs")

vars = [
    "phi", 
    "sumE", 
    # "maxE", 
    # "minE",
    # "avgE", 
    "r", 
    "singleE",
    "dy", 
    "dz",
    "sameSide",
    "lPint",
    "lPext",
    "trackLength1",
    "trackLength2",
    # "caloTime1",
    # "caloTime2",
    ]

bins = (
    phi = (0,180),
    sumE = (0, 3500),
    # minE = (0, 3500),
    # maxE = (0, 3500),
    # avgE = (0, 3500),
    r = (0, 200),
    singleE = (0, 3500),
    dy = (0, 200),
    dz = (0, 200),
    sameSide = (0, 2),
    lPint = (0, 50),
    lPext = (0, 50),
    trackLength1 = (0, 3000),
    trackLength2 = (0, 3000),
    # caloTime1 = (0, 100),
    # caloTime2 = (0, 100),
)

processes = load_ndim_processes(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), bins, vars)

data_dir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/neutrons_jan_2026/")

include(joinpath(data_dir, "read_neutrons_ND.jl"))
neutron_processes = load_neutron_processes_ND(data_dir, vars, bins)
for p in neutron_processes
    set_activity!(p, p.activity / 3)
    println("Neutron process: ", p.isotopeName, " nTotalSim: ", p.nTotalSim, " activity: ", p.activity)
end


signal = get_process(signal_name, processes) |> first

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
if radon_tag == 1
    set_activity!(background[7], 150/1e6) # radon to 150 uBq/kg
elseif radon_tag == 2
    set_activity!(background[7], 2/1e3) # radon to 2 mBq/kg
elseif radon_tag == 3
    set_activity!(background[7], 0.6/1e3) # radon to 0.6 mBq/kg
end

α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")

function filter_process!(process::DataProcessND, var_name::Symbol, range::Tuple{Real, Real})
    process.data = filter(x -> getproperty(x, var_name) >= range[1] && getproperty(x, var_name) <= range[2], process.data)
    @info "filtered process $(process.isotopeName) on $var_name in range $range, new nEvents: $(length(process.data))"
    return process;
end

side_range = ifelse(side == "sameSide", (1., 1.), (0.,0.))
filter_process!(signal, :sameSide, side_range);
map( x-> filter_process!(x, :sameSide, side_range), background);




prob(x) = - get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="formula")


function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


searchRange = make_stepRange(signal)


lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float

############### Optimization
using Metaheuristics
options = Options(;
    # x_tol = 1.0,
    # f_tol = 1e-5,
    # f_tol_rel = 1e-5,
    # f_tol_abs = 1e-5,
    time_limit = 60*60*22.0,
    # parallel_evaluation = true,
    verbose = true,
    # iterations = 15,
    # store_convergence = true
)

bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

function f_parallel(X)
    fitness = zeros(size(X,1))
    Threads.@threads for i in 1:size(X,1)
        fitness[i] = prob(X[i,:])
    end
    fitness
end

x0 = float.([
    5, 180, # phi
    2600, 3000, # sumE
    0, 3500, # singleE
    0,150, # r
    0, 150, # dy
    0, 150, # dz
    0, 2, # sameSide
    0, 3, # lPint
    1, 100, # lPext
    0, 2500, # trackLength1
    0, 2500, # trackLength2
    # 0, 20, # caloTime1
    # 0, 20, # caloTime2
])


algo = PSO(;options)
set_user_solutions!(algo, x0, prob)


# result = Metaheuristics.optimize(f_parallel, bounds, algo)
# result = Metaheuristics.optimize(prob, bounds, PSO(;options))
result = Metaheuristics.optimize(prob, bounds, algo)
@show minimum(result)
@show res=  minimizer(result)

function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => (best[i], best[i+1]) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
    best_roi = NamedTuple(
        k => (best[i], best[i+1]) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

best = get_best_ROI_ND(res, signal)
best_sens = get_sensitivityND(SNparams, α, vcat(signal, background), best; approximate="table")

import Base.print
function print(best_roi::SensitivityEstimateND)
    println("tHalf: $(best_roi.tHalf)")
    println("signalEff: $(best_roi.signalEff)")
    println("bkg count: $(best_roi.bkgCounts)")
    println("best ROI:")
    for (k,v) in pairs(best_roi.roi)
        println("  $k : $(v)")
    end
end

print(best_sens)
@info print("Best sensitivity found:", best_sens)

best = get_best_ROI_ND(res, signal)

dff = DataFrame(
    tHalf = best_sens.tHalf,
    eff = best_sens.signalEff,
    best_roi = best_sens.roi
)

open(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/nd_cc_runs/results/20260123_results/$(now)_$(side)_$(signal_name)_$(rand(1:1000000))_radon_config_$(radon_tag).csv"), "w") do io
    CSV.write(io, dff; header=true, delimiter=',')
end

