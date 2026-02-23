using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

# push!(LOAD_PATH, srcdir())
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
# include(srcdir("params/Params.jl"))


## Helper functions
include(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/helper_functions.jl"))

vars = [
    "phi", 
    "sumE", 
    # "maxE", 
    # "minE",
    # "avgE", 
    "r", 
    # "singleE",
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
    # singleE = (0, 3500),
    dy = (0, 200),
    dz = (0, 200),
    sameSide = (0, 1.0),
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
background = vcat(background, neutron_processes)

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)
# set_activity!(background[7], 0.15/1e3) # radon to 0.15mBq/kg
# set_activity!(background[7], 2/1e3) # radon to 2mBq/kg
radon_tag = 2
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

side = "sameSide" # "oppositeSide
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
    time_limit = 10*60.0, #60*60*12.0,
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
    500, 3000, # sumE
    # 0, 3500, # singleE
    0,150, # r
    0, 150, # dy
    0, 150, # dz
    0,2, # sameSide
    0, 1, # lPint
    1, 100, # lPext
    0, 2500, # trackLength1
    0, 2500, # trackLength2
    # 0, 20, # caloTime1
    # 0, 20, # caloTime2
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
@info "Best sensitivity found: $(best_sens.tHalf) years"

best = get_best_ROI_ND(res, signal)
h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external = get_background_counts_hists(background, best, 0:100:3500, :sumE)

@show integral.([ h_bb,h_internal, h_radon, h_Bi210, h_detector, h_external])

let 
    my_roi = (
        phi= (0, 180),
        sumE= (2700, 3000),
        # singleE= (0, 2600),
        r= (0, 200),
        dy= (0, 70),
        dz= (0, 80),
        sameSide = (0.0, 2.0),
        lPint= (0, 3.), # 5.0),
        lPext= (1., 50),
        trackLength1= (0, 3000),
        trackLength2= (0, 3000),
    )
    my_sens = get_sensitivityND(SNparams, α, vcat(signal, background), my_roi; approximate="table")
    print(my_sens)
end



dff = DataFrame(
    tHalf = best_sens.tHalf,
    eff = best_sens.signalEff,
    best_roi = best_sens.roi
)

# open("/pbs/home/m/mpetro/sps_mpetro/Projects/PhD/SNSensitivityEstimate/scripts/ND/results/Jan26/$(signal_name)_$(rand(1:1000000)).csv", "w") do io
#     CSV.write(io, dff; header=true, delimiter=',')
# end


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
# colors = ColorSchemes.roma100






begin
    var = :phi
    binning = 0:1:180
    roi = (
        phi= (0, 180),
        sumE= (0, 3500),
        r= (0, 200),
        dy= (0, 200),
        dz= (0, 200),
        sameSide = (0.0, 2.0),
        lPint= (0.0, 50),
        lPext= (0, 50.0),
        trackLength1= (0, 5000),
        trackLength2= (0, 5000),
    )

    h_internal, h_radon, h_Bi10, h_detector, h_external = get_background_counts_hists(background, roi, binning, var)

    f = Figure(size = (1000, 1400), fontsize = 32)
    a1 = Axis(f[1,1], title = "signal", tellwidth = false, width = 800) 
    a2 = Axis(f[2,1], title = "internal", tellwidth = false, width = 800) 
    a3 = Axis(f[3,1], title = "radon", tellwidth = false, width = 800, ) 
    a4 = Axis(f[4,1], title = "Bi210", tellwidth = false, width = 800, ) 
    a5 = Axis(f[5,1], title = "detector", tellwidth = false, width = 800) 
    a6 = Axis(f[6,1], title = "external", tellwidth = false, width = 800, ) 
    ps = []
    Label(f[0,1], "distributions for variable: $var", fontsize = 36)

    hidexdecorations!.([a1, a2, a3, a4, a5])

    # d_signal = log10.(abs.((getproperty.(signal.data, :caloTime2) .- getproperty.(signal.data, :caloTime1))) ./ (getproperty.(signal.data, :sumE)))
    d_signal = getproperty.(signal.data, var) 
    # p = stephist!(a, d_signal, bins=binning, label = "signal", color =colors[1], offset = 0, scale_to = 0.7)
    p = hist!(a1, d_signal, bins=binning,  color =colors[1])
    
    push!(ps, p)
    
    labels = ["signal", "internal", "radon", "Bi210", "detector", "external"]
    
    push!(
        ps,
        stephist!(a2, h_internal; color = colors[3], linewidth = 4)
    )
    push!(
        ps,
        stephist!(a3, h_radon; color = colors[5], linewidth = 4)
    )
    push!(
        ps,
        stephist!(a4, h_Bi10; color = colors[7], linewidth = 4)
    )
    push!(
        ps,
        stephist!(a5, h_detector; color = colors[9], linewidth = 4)
    )

    push!(
        ps,
        stephist!(a6, h_external; color = colors[11], linewidth = 4)
    )

    vlines!.([a1, a2, a3, a4, a5, a6], [2.0],color = :black, linestyle = :solid, label = "data-cut line", linewidth = 5)
    axislegend(a1, position = :rt, padding = (10,10,10,10))
    # Legend(f[1:end,2], ps, labels, tellheight = true)
    # save(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/figs/ND_backgrounds_$(var).png"), f, px_per_unit = 2)

    f
end


colors = ["#041E42","#BE4D00","#951272","#006630","#005C8A","#FFB948","#605643","#302D23"]

roi
set_activity!(background[7], 2/1e3) # radon to 2mBq/kg

begin
    begin
        binning = 300:100:3500
        bw = step(binning)
        var = :sumE
        fill_area = true
        stacked = true
        roi = (
            phi= (0, 180),
            sumE= (300, 3500),
            r= (0, 100),
            dy= (0, 100),
            dz= (0, 100),
            lPint= (0.0, prob_to_lprob(0.04)),
            lPext= (prob_to_lprob(0.01), 50.0),
            trackLength1= (0, 2500),
            trackLength2= (0, 2500),
        )
        scale = :identity
        colors = colors
    end

    # h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external = get_background_counts_hists(background, roi, binning, var)

    @show integral.([ h_bb,h_internal, h_radon, h_Bi210, h_detector, h_external])

    f = Figure(size = (1900, 1100), fontsize = 36, figure_padding = 16)

    sum_hists = [h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external]
   
    legend_elements = [PolyElement(polycolor = colors[i]) for i in 1:length(sum_hists)] 
    legend_labels = [
        L"2\nu\beta\beta",
        "Internal contamination",
        "Radon",
        "Bi-210 wire surface",
        "Detector contamination",
        "External",
    ]

    legend_master = Legend(
        f, legend_elements, legend_labels,
        patchsize = (45, 35),
        rowgap = 20,
        tellheight = true,
        padding = (10, 10, 10, 10),
        nbanks = 1,
    )

    legend_master.nbanks = 1
    legend_master.tellheight = false
    legend_master.tellwidth = true

    title = "Stacked histogram\nSuperNEMO expected background for 17.5kg.yr exposure\nvariable: $var"


    ax = Axis(
        f[1, 1];
        xlabel = L"Summed 2-electron energy (keV)$$",
        # xlabel = L"$\varphi$ (°)$$",
        ylabel = L"Counts / (\mathrm{yr} \cdot %$bw \mathrm{keV})$$",
        title  = title,
        yscale = scale == :log10 ? log10 : identity,
        yminorticks=IntervalsBetween(10), 
        yminorticksvisible = true, 
        yminorticksize = 4
    )
    hist!(ax, sum(sum_hists), color = (:black, 0.15), label = "Total background")
    if !fill_area
        axislegend(ax, position = :rt, padding = (20,20,20,20))
    end

    if stacked
        stacked_hist!(ax, sum_hists; labels = nothing, colors = colors, fill_area = fill_area)
    else
        all_individual_hist!(ax, sum_hists; labels = nothing, colors = colors, fill_area = fill_area)
    end

    min_cts = max(1e-3,minimum(bincounts(sum(sum_hists))))
    max_cts = find_max_bincounts(sum(sum_hists))

    if scale == :log10
        ylims!(ax, 10^-3, 1e2 * max_cts)
        decades = -3:ceil(log10(1e2*max_cts)) .|> Int
        ticks   = 10.0 .^ decades 
        labels_  = [L"10^{%$d}" for d in decades]
        ax.yticks = ( ticks,labels_ )
        

    else
        ylims!(ax, 0, 1.2*max_cts)  
    end

    xlims!(ax, first(binning), last(binning))

    f[1, 2] = legend_master

    save(joinpath(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/figs"),
        "bad_radon_mode_$(var)_background_categories_$(scale)_fill_$(fill_area)_stacked_$(stacked).png"), f)

    f
end






