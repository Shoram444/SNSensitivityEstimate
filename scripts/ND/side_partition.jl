using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie, Surrogates, Revise
Revise.track(SensitivityModule)


# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(scriptsdir("Params.jl"))

vars = [
    "phi", 
    # "sumE", 
    # "maxE", 
    # "minE", 
    # "r", 
    # "dy", 
    # "dz",
    "sameSide"
    ]

bins = (
    phi = (0,180),
    # sumE = (0, 3500),
    # maxE = (0, 3500),
    # minE = (0, 3500),
    # r = (0, 100),
    # dy = (-100, 100),
    # dz = (-100, 100),
    sameSide = (0, 1)
)

processes = load_ndim_processes("fal5_TKrec", bins, vars)

"""
    Function that partitions the `DataProcessND` object to n objects based on the key (i.e. sameSide) and the keyValues (i.e. 0, 1).
    The function creates a new `DataProcessND` object for each partition and appends it to the partitions array.
    The function also checks if the key is present in the varNames of the `DataProcessND` object. If not, it throws an error.
    The function returns an array of `DataProcessND` objects.
    # Arguments
    - `processes::Vector{<:DataProcessND}`: Vector of `DataProcessND` objects to be partitioned.
    - `key::String`: The key to partition the `DataProcessND` object. Default is "sameSide".
    - `keyValues::Vector{<:Real}`: The values of the key to partition the `DataProcessND` object. Default is [0, 1].
    # Returns
    - `partitions::Vector{<:DataProcessND}`: Vector of `DataProcessND` objects that are partitioned based on the key and keyValues.
    # Example
    ```julia
    processes = load_ndim_processes("fal5_TKrec", bins, vars)
    partitions = make_partitions(processes, "sameSide", [0, 1])
    ```
"""
function make_partitions(processes, key = "sameSide", keyValues = [0, 1])
    findfirst(x -> x == key, processes[1].varNames) === nothing && @error "sameSide not found in varNames"

    partitions = Partition[]

    for p in processes
        for k in keyValues
            new_p = deepcopy(p)
            new_p.data = filter(x -> getproperty(x, Symbol(key)) == k, p.data)
            push!(partitions, Partition(new_p, key, k))
        end
        println("partitions for $(p.isotopeName) created!")
    end
    return partitions
end

mutable struct Partition
    process::DataProcessND
    key::String
    keyValue::Real
end

partitions = make_partitions(processes)

# signal = get_process("bb0nu_foil_bulk", processes)
signal = get_process("RH037_foil_bulk", partitions)
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
bkgs = [
    get_process("bb_foil_bulk", partitions),
    get_process("Bi214_foil_bulk", partitions),
    get_process("Bi214_wire_surface", partitions),
    get_process("Tl208_foil_bulk", partitions),
    get_process("K40_foil_bulk", partitions),
    get_process("Pa234m_foil_bulk", partitions),
]

background = vcat(bkgs...)

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)
set_signal!(background[2], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!.( signal, 1e8 )
# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!.( background[1:2], 0.99e8 )
set_nTotalSim!.( background[3:4], 1e8 )
set_nTotalSim!.( background[5:6], 1e8 )
set_nTotalSim!.( background[7:8], 1e8 )
set_nTotalSim!.( background[9:10], 1e8 )
set_nTotalSim!.( background[11:12], 1e8 )

α= 1.64485362695147


using FHist, CairoMakie
let 
    h1 = Hist1D(getproperty.(signal[1].data, :phi); binedges=0:180)
    h2 = Hist1D(getproperty.(signal[2].data, :phi); binedges=0:180)

    f = Figure()
    ax = Axis(f[1, 1], title="phi", xlabel="phi", ylabel="counts", )
    plot!(ax, h1, label = "partition 1")
    plot!(ax, h2, label = "partition 2")
    f
    
end

α= 1.64485362695147

prob(x) = - SensitivityModule.get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="formula")


function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


searchRange = make_stepRange(signal[1])

lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float

############### Optimization
using Metaheuristics
options = Options(;
    x_tol = 1.0,
    f_tol = 1e-1,
    f_tol_rel = 1e-1,
    f_tol_abs = 1e-1,
    time_limit = Inf,
    parallel_evaluation = true,
    verbose = true,
    iterations = 2,
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
