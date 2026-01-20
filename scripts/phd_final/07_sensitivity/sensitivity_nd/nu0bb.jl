using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

# push!(LOAD_PATH, srcdir())
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
# include(srcdir("params/Params.jl"))

vars = [
    "phi", 
    "sumE", 
    "maxE", 
    "minE",
    "avgE", 
    "r", 
    "singleE",
    "dy", 
    "dz",
    "sameSide",
    "lPint",
    "lPext",
    "trackLength1",
    "trackLength2",
    "caloTime1",
    "caloTime2",
    ]

bins = (
    phi = (0,180),
    sumE = (2600, 3500),
    minE = (0, 3500),
    maxE = (0, 3500),
    avgE = (0, 3500),
    r = (0, 200),
    singleE = (0, 3500),
    dy = (0, 200),
    dz = (0, 200),
    sameSide = (0, 1),
    lPint = (0, 100),
    lPext = (0, 100),
    trackLength1 = (0, 5000),
    trackLength2 = (0, 5000),
    caloTime1 = (0, 1000),
    caloTime2 = (0, 1000),
)

processes = load_ndim_processes(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), bins, vars)

# signal_name = "%SIGNAL"
signal_name = "bb0nu_foil_bulk"
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

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)


α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


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
    time_limit = 1*20.0, #60*60*12.0,
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
    0, 180, # phi
    2700, 3100, # sumE
    0, 3500, # singleE
    0, 100, # dy
    0, 100, # dz
    0, 3, # lPint
    1, 100, # lPext
    0, 2000, # trackLength1
    0, 2000, # trackLength2
])


algo = PSO(;options)
# set_user_solutions!(algo, reshape(x0, 1, :), prob)
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

let 
    my_roi = (
        phi= (5, 180),
        sumE= (2700, 3100),
        singleE= (0, 3500),
        dy= (0, 100),
        dz= (0, 100),
        lPint= (0, 3.0), # 5.0),
        lPext= (1.0, 100),
        trackLength1= (0, 2000),
        trackLength2= (0, 2000),
    )
    my_sens = get_sensitivityND(SNparams, α, vcat(signal, background), my_roi; approximate="table")
    print(my_sens)
end


dff = DataFrame(
    tHalf = best_sens.tHalf,
    eff = best_sens.signalEff,
    best_roi = best_sens.roi
)

open("/pbs/home/m/mpetro/sps_mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND/results/May25_tof_r/$(signal_name)_$(rand(1:1000000)).csv", "w") do io
    CSV.write(io, dff; header=true, delimiter=',')
end


using FHist


# pint = getproperty.(signal.data, :lPint) 
# # replace!(pint, -Inf32 => 110)
# pext = getproperty.(signal.data, :lPext) 
# # replace!(pext, -Inf32 => 110)

# pintb = getproperty.(background[end].data, :lPint) 
# # replace!(pintb, -Inf32 => 110)
# pextb = getproperty.(background[end].data, :lPext) 
# # replace!(pextb, -Inf32 => 110)

# begin
#     name = "gamma"
#     f = Figure(size = (1800,800), fontsize = 25)
#     a = Axis(f[1,1], xlabel = "Pint", ylabel = "Pext", title = "signal: -log(P_tof)")
#     a2 = Axis(f[1,3], xlabel = "Pint", ylabel = "Pext", title = "$name: -log(P_tof)")

#     colorscale = cgrad(:plasma, 1, scale = :log10)
#     colorscaleb = cgrad(:plasma, 1, scale = :log10)

#     h2 = Hist2D((abs.(pint), abs.(pext)); binedges = (0:2:110, 0:2:110)) 
#     h2b = Hist2D((abs.(pintb), abs.(pextb)); binedges = (0:2:110, 0:2:110)) 
#     p = plot!(a, h2, colormap = colorscale)

#     p2 = plot!(a2, h2b, colormap = colorscaleb)
#     Colorbar(f[1,2], p, label = "normalized counts")
#     Colorbar(f[1,4], p2, label = "normalized counts")
#     save( "Pint_Pext_signal_$name.png", f, px_per_unit = 2)
#     f
# end

using ColorSchemes
colors = ColorSchemes.tab20

begin
    var = "log10((:caloTime2 - :caloTime1) / :sumE)"
    binning = -5:0.01:5
    f = Figure(size = (1200, 900), fontsize = 25)
    a = Axis(f[1,1], xlabel = string(var), ylabel = "Normalized counts", title = "ND Backgrounds and Signal - variable: $var") 
    ps = []

    d_signal = log10.(abs.((getproperty.(signal.data, :caloTime2) .- getproperty.(signal.data, :caloTime1))) ./ (getproperty.(signal.data, :sumE)))
    # d_signal = getproperty.(signal.data, var) .- 
    p = hist!(a, d_signal, bins=binning, label = "signal", color =colors[1], offset = length(background)+1, scale_to = 0.7)
    
    push!(ps, p)
    
    labels = ["signal"]
    labels = append!(labels, [b.isotopeName for b in background] )
    
    for (i, j) in zip(length(background):-1:1,1:length(background))
        d = log10.(abs.((getproperty.(background[j].data, :caloTime2) .- getproperty.(background[j].data, :caloTime1))) ./ (getproperty.(background[j].data, :sumE)))
        # d = getproperty.(background[j].data, var)
        # h = normalize(Hist1D(d; binedges=binning))
        p = hist!(a, d, bins=binning, color = colors[j+1], label = background[j].isotopeName,offset = i, scale_to = 0.7)
        push!(ps, p)
    end
    # stephist!(a, normalize(Hist1D(getproperty.(background[1].data, :singleE); binedges=0:50:3500)), label = "2nubb", linewidth = 3)
    # stephist!(a, normalize(Hist1D(getproperty.(background[2].data, :singleE); binedges=0:50:3500)), label = "Bi214", linewidth = 3)
    # stephist!(a, normalize(Hist1D(getproperty.(background[5].data, :singleE); binedges=0:50:3500)), label = "Tl208", linewidth = 3)
    # stephist!(a, normalize(Hist1D(getproperty.(background[10].data, :singleE); binedges=0:50:3500)), label = "K40 PMT", linewidth = 3)
    # axislegend(a)
    Legend(f[1,2], ps, labels)

    f
end