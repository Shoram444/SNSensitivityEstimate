using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie 
using Revise
Revise.track(SensitivityModule)

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(srcdir("params/Params.jl"))

vars = [
    "phi", 
    "sumE", 
    # "maxE", 
    # "minE", 
    "r", 
    # "dy", 
    # "dz",
    # "sameSide"
    ]

bins = (
    phi = (0,180),
    sumE = (0, 3500),
    # maxE = (0, 3500),
    # minE = (0, 3500),
    r = (0, 100),
    # dy = (-100, 100),
    # dz = (-100, 100),
    # sameSide = (0, 1)
)

processes = load_ndim_processes("fal5_TKrec", bins, vars)

# signal = get_process("bb0nuM1_foil_bulk", processes) |> first
signal = get_process("RH037_foil_bulk", processes) |> first
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", processes) |> first,
    get_process("Bi214_foil_bulk", processes) |> first,
    get_process("Bi214_wire_surface", processes) |> first,
    get_process("Tl208_foil_bulk", processes) |> first,
    get_process("K40_foil_bulk", processes) |> first,
    get_process("Pa234m_foil_bulk", processes) |> first,
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!( signal, 1e8 )
# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 1e8 )
set_nTotalSim!( background[2], 1e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 1e8 )
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )

α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


prob(x) = - SensitivityModule.get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="table")


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
    x_tol = 1.0,
    f_tol = 1e-1,
    f_tol_rel = 1e-1,
    f_tol_abs = 1e-1,
    time_limit = 60*60*1.0,
    parallel_evaluation = true,
    verbose = true,
    iterations = 15,
    store_convergence = true
)

bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

function f_parallel(X)
    fitness = zeros(size(X,1))
    Threads.@threads for i in 1:size(X,1)
        fitness[i] = prob(X[i,:])
    end
    fitness
end

result = Metaheuristics.optimize(f_parallel, bounds, ECA(;options))
# result = Metaheuristics.optimize(prob, bounds, SA(;options))
@show minimum(result)
@show res=  minimizer(result)

function get_best_ROI_ND(res, process)
    best = best_candidate(res)
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

best = get_best_ROI_ND(res, signal)
get_sensitivityND(SNparams, α, vcat(signal, background), best; approximate="table")

res2 = float.([10,175,2710,3350,0,50])
best2 = get_best_ROI_ND(res2, signal)
get_sensitivityND(SNparams, α, vcat(signal, background), best2; approximate="table", add_mock_bkg=1.43)


begin
    res2 = float.([30,175, 0, 3200, 0, 50])
    best2 = get_best_ROI_ND(res2, signal)
    get_sensitivityND(SNparams, α, vcat(signal, background), best2; approximate="table", add_mock_bkg=1.2)
end


# f_calls, best_f_value = convergence(result)
# plot(f_calls, best_f_value,)

