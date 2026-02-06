using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist

println("loaded pkgs")


## Helper functions
include(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/helper_functions.jl"))

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
    # "sameSide",
    "Pint",
    "Pext",
    "lPint",
    "lPext",
    "trackLength1",
    "trackLength2",
    "caloTime1",
    "caloTime2",
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
    # sameSide = (0, 1),
    Pint = (0, 1),
    Pext = (0, 1),
    lPint = (0, 50),
    lPext = (0, 50),
    trackLength1 = (0, 3000),
    trackLength2 = (0, 3000),
    caloTime1 = (0, 100),
    caloTime2 = (0, 100),
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
radon_tag = 1
if radon_tag == 1
    set_activity!(background[7], 150/1e6) # radon to 150 uBq/kg
elseif radon_tag == 2
    set_activity!(background[7], 2/1e3) # radon to 2 mBq/kg
elseif radon_tag == 3
    set_activity!(background[7], 0.6/1e3) # radon to 0.6 mBq/kg
end

α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


prob(x) = - get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="formula")


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
    5, 180, # phi
    2700, 3000, # sumE
    # 0, 3500, # singleE
    0,150, # r
    0, 150, # dy
    0, 150, # dz
    0, 3, # lPint
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
# h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external = get_background_counts_hists(background, best, 0:100:3500, :sumE)

# @show integral.([ h_bb,h_internal, h_radon, h_Bi210, h_detector, h_external])

let 
    var = :sumE
    binning = 0:100:3500
    my_roi = (
        phi = (0.0, 170.0), 
        sumE = (400, 2700), 
        r = (0.0, 200.0), 
        singleE = (0.0, 2600), 
        dy = (0.0, 200.0), 
        dz = (0.0, 200.0), 
        lPint = (0.0, 3), 
        lPext = (1.0, 50.0), 
        trackLength1 = (0.0, 2200.0), 
        trackLength2 = (0.0, 2200.0)
    )
    my_sens = get_sensitivityND(SNparams, α, vcat(signal, background), my_roi; approximate="table")
    print(my_sens)

    # h_internal, h_radon, h_Bi10, h_detector, h_external = get_background_counts_hists(background, my_roi, binning, var);
    # @show integral.([ h_internal, h_radon, h_Bi10, h_detector, h_external])
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


d_signal = (abs.((getproperty.(signal.data, :caloTime2) .- getproperty.(signal.data, :caloTime1))) ./ sqrt.((getproperty.(signal.data, :sumE))))


# plots of different variable distributions for each process
begin
    var = :lPint
    binning = 0:0.1:50
    roi = (
        phi= (0, 180),
        sumE= (300, 3500),
        r= (0, 100),
        singleE = (0, 3500),
        dy= (0, 100),
        dz= (0, 100),
        Pint = (0, 1),
        Pext = (0, 1),
        lPint= (0.0, 50), #prob_to_lprob(0.04)),
        lPext= (0, 50), #(prob_to_lprob(0.01), 50.0),
        trackLength1= (0, 2500),
        trackLength2= (0, 2500),
        caloTime1 = (0, 100),
        caloTime2 = (0, 100),
    )

    h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external = get_background_counts_hists(background, roi, binning, var)

    f = Figure(size = (900, 1400), fontsize = 16)
    a1 = Axis(f[1,1], title = L"$$signal", ylabel = "counts (a.u.)",tellwidth = false, width = 650) 
    a2 = Axis(f[2,1], title = L"2\nu\beta\beta", ylabel = "counts / (17.5 kgy)", tellwidth = false, width = 650)
    a3 = Axis(f[3,1], title = L"$$internal", ylabel = "counts / (17.5 kgy) ", tellwidth = false, width = 650) 
    a4 = Axis(f[4,1], title = L"$$radon", ylabel = "counts / (17.5 kgy) ", tellwidth = false, width = 650, ) 
    a5 = Axis(f[5,1], title = L"${}^{210}$Bi", ylabel = "counts / (17.5 kgy) ", tellwidth = false, width = 650, ) 
    a6 = Axis(f[6,1], title = L"$$detector", ylabel = "counts / (17.5 kgy) ", tellwidth = false, width = 650) 
    a7 = Axis(f[7,1], title = L"$$external", ylabel = "counts / (17.5 kgy) "
        # , xlabel = L"$\varphi$ (°)", 
        # , xlabel = L"$E_{sum}$ (keV)", 
        # , xlabel = L"$r$ (mm)", 
        # , xlabel = L"$log(P_{int})$", 
        , xlabel = L"$log(P_{ext})$", 
        tellwidth = false, width = 650, ) 
    ps = []
    Label(f[0,1], "distributions for variable: $var", fontsize = 26)
    # Label(f[0,1], "distributions for variable: caloTime2 - caloTime1 / sqrt(sumE)", fontsize = 36)

    hidexdecorations!.([a1, a2, a3, a4, a5, a6])

    # d_signal = (abs.((getproperty.(signal.data, :caloTime2) .- getproperty.(signal.data, :caloTime1))) ./ sqrt.((getproperty.(signal.data, :sumE))))
    d_signal = getproperty.(signal.data, var) 
    # p = stephist!(a, d_signal, bins=binning, label = "signal", color =colors[1], offset = 0, scale_to = 0.7)
    p = hist!(a1, d_signal, bins=binning,  color =colors[1])


    push!(ps, p)

    # d2 = (abs.((getproperty.(background[1].data, :caloTime2) .- getproperty.(background[1].data, :caloTime1))) ./ sqrt.((getproperty.(background[1].data, :sumE))))
    # p2 = hist!(a2, d2, bins=binning, color =colors[1+1])

    labels_ = ["signal", L"2\nu\beta\beta","internal", "radon", L"${}^{210}$Bi", "detector", "external"]
    
    push!(
        ps,
        # p2
        stephist!(a2, h_bb; color = colors[2], linewidth = 4)
    )
    push!(
        ps,
        # p2
        stephist!(a3, h_internal; color = colors[3], linewidth = 4)
    )
    # d3 = (abs.((getproperty.(background[7].data, :caloTime2) .- getproperty.(background[7].data, :caloTime1))) ./ sqrt.((getproperty.(background[7].data, :sumE))))
    # p3 = hist!(a3, d3, bins=binning, color =colors[7+1])

    push!(
        ps,
        # p3
        stephist!(a4, h_radon; color = colors[5], linewidth = 4)
    )
    # d4 = (abs.((getproperty.(background[6].data, :caloTime2) .- getproperty.(background[6].data, :caloTime1))) ./ sqrt.((getproperty.(background[6].data, :sumE))))
    # p4 = hist!(a4, d4, bins=binning, color =colors[5+1])
    push!(
        ps,
        # p4
        stephist!(a5, h_Bi210; color = colors[7], linewidth = 4)
    )
    push!(
        ps,
        stephist!(a6, h_detector; color = colors[9], linewidth = 4)
    )
    # d_external = []
    # for b in 19:length(background)
    #     append!(d_external, (abs.((getproperty.(background[b].data, :caloTime2) .- getproperty.(background[b].data, :caloTime1))) ./ sqrt.((getproperty.(background[b].data, :sumE)))))
    # end
    # p6 = hist!(a6, d_external, bins=binning, color =colors[9+1])
    push!(
        ps,
        # p6,
        stephist!(a7, h_external; color = colors[11], linewidth = 4)
    )


    # vlines!.([a1, a2, a3, a4, a5, a6], [2.0],color = :black, linestyle = :solid, label = "data-cut line", linewidth = 5)
    # axislegend(a1, position = :rt, padding = (10,10,10,10))
    # Legend(f[1:end,2], ps, labels_, tellheight = true)
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/figs/ND_backgrounds_$(var).png"), f, px_per_unit = 2)

    f
end



colors = ["#041E42","#BE4D00","#951272","#006630","#005C8A","#FFB948","#605643","#302D23"]

roi
set_activity!(background[7], 2/1e3) # radon to 2mBq/kg

# total background figure
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
            singleE = (0, 3500),
            dy= (0, 100),
            dz= (0, 100),
            Pint = (0, 1),
            Pext = (0, 1),
            lPint= (0.0, 50), #prob_to_lprob(0.04)),
            lPext= (0, 50), #(prob_to_lprob(0.01), 50.0),
            trackLength1= (0, 2500),
            trackLength2= (0, 2500),
            caloTime1 = (0, 100),
            caloTime2 = (0, 100),
        )
        scale = :identity
        colors = colors
    end

    h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external = get_background_counts_hists(background, roi, binning, var)

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





# 2D figure for multidim roi
# begin
    var1 = :sumE
    # var2 = :lPint
    var2 = :phi
    binning1 = 0:10:3500
    # binning2 = 0:0.5:35.0
    binning2 = 0:0.1:180
    # roi = (
    #     phi= (0, 180),
    #     sumE= (0, 3500),
    #     r= (0, 200),
    #     singleE = (0, 2650),
    #     dy= (0, 100),
    #     dz= (0, 100),
    #     Pint = (0, 1),
    #     Pext = (0, 1),
    #     lPint= (0.0, 4.0), #prob_to_lprob(0.04)),
    #     lPext= (2.0, 50), #(prob_to_lprob(0.01), 50.0),
    #     trackLength1= (0, 3000),
    #     trackLength2= (0, 3000),
    #     caloTime1 = (0, 100),
    #     caloTime2 = (0, 100),
    # )

    roi = (
        phi = (0, 180), 
        sumE = (0, 3500), 
        dy = (0, 135), 
        dz = (0, 140), 
        lPint = (0.0, 4.5), 
        lPext = (2.4, 50.0)
    )

    # h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external = get_background_counts_hists(background, roi, binning1, var1)

    # sum(integral.([ restrict(b, 2700, 300) for b in [h_bb,h_internal, h_radon, h_Bi210, h_detector, h_external]]))

    # tot_bkg = sum(integral.([ h_bb,h_internal, h_radon, h_Bi210, h_detector, h_external]))
    function apply_roi_to_process!(p, roi)
        data = p.data
        varNames = keys(roi)
        if roi !== nothing
            for (i,n) in enumerate(varNames)
                data = filter(x -> getproperty(x, Symbol(n)) > roi[i][1] && getproperty(x, Symbol(n)) < roi[i][2], data)
            end
        end
        p.data = data
        return p
    end

    for p in vcat(signal, background)
        apply_roi_to_process!(p, roi)
    end

    tot_signal_eff = length(signal.data) / signal.nTotalSim

    d1_s = getproperty.(signal.data, var1)
    d2_s = getproperty.(signal.data, var2)
    h_signal = Hist2D(( d1_s, d2_s ); binedges = (binning1, binning2)) |> normalize
    h_signal = h_signal * tot_signal_eff

    h_background = Hist2D(; binedges = (binning1, binning2)) 

    # figure out wtf is wrong with externals! 
    for b in background
        d1_b = getproperty.(b.data, var1)
        d2_b = getproperty.(b.data, var2)
        h_b = Hist2D(( d1_b, d2_b ); binedges = (binning1, binning2)) |> normalize

        b_cts = integral(get_roi_bkg_counts_hist(b, roi, binning1, var1))
        isotope = b.isotopeName
        println("Background process: $isotope total counts in ROI: $b_cts")
        h_b = h_b * b_cts
        b_cts > 0 && (h_background += h_b)
    end


    # d1_b = vcat([getproperty.(b.data, var1) for b in background]...)
    # d2_b = vcat([getproperty.(b.data, var2) for b in background]...)
    # h_background = Hist2D(( d1_b, d2_b ); binedges = (binning1, binning2)) |> normalize
    # h_background = h_background * tot_bkg

    f = Figure(size = (2400,800), fontsize = 38)
    # a1 = Axis(f[1,1], xlabel = L"Summed 2-electron energy (keV) $$", ylabel = L"$\varphi$ (°) $$", title = L"signal: $0\nu\beta\beta$", )
    a1 = Axis(f[1,1], xlabel = L"Summed 2-electron energy (keV) $$", ylabel = L"$-log(P_\mathrm{int})$", title = L"signal: $0\nu\beta\beta$", )
    # a2 = Axis(f[1,3], xlabel = L"Summed 2-electron energy (keV) $$", ylabel = L"$\varphi$ (°) $$", title = L"total background $$", )
    a2 = Axis(f[1,3], xlabel = L"Summed 2-electron energy (keV) $$", ylabel = L"$-log(P_\mathrm{int})$", title = L"total background $$", )

    # colorscale = cgrad(:viridis, 1, scale = :log10)
    colormap = [Makie.to_color(:transparent); Makie.to_colormap(:viridis)]
    p1= plot!(a1, h_signal, colorscale = log10, colormap = colormap)
    p2 =plot!(a2, h_background, colorscale = log10 , colormap = colormap)
    Colorbar(f[1,2], p1, label = "signal efficiency")
    Colorbar(f[1,4], p2, label = "expected background")

    poly!(
		a1,
		[2700,2700,3000,3000],
		[0,4,4,0],
		color = (:black,0),
		strokewidth = 5,
		strokecolor = :red,
	 )
	poly!(
		a2,
		[2700,2700,3000,3000],
		[0,4,4,0],
		color = (:black,0),
		strokewidth = 5,
		strokecolor = :red,
	 )
    colgap!(f.layout, 2, Relative(0.05))
    # save( joinpath(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/figs"),
    #     "ND_2Dhist_signal_background_$(var1)_vs_$(var2).png"), f, px_per_unit = 2)
    f


# end

# extract from the signal histogram a list of midpoints for binedge1, binedge2 and corresponding counts
binedges1 = midpoints(binedges(h_signal)[1]) |> collect
binedges2 = midpoints(binedges(h_signal)[2]) |> collect

# make combinations of all binedges by taking the midpoints so we have x,y coordinate for each bin where x, y are the midpoints
bin1_bin2 = [ (bin_mid_1, bin_mid_2) for bin_mid_1 in binedges1 for bin_mid_2 in binedges2 ]


# get counts for each bin
counts = [lookup(h_signal, bin1, bin2) for (bin1, bin2) in bin1_bin2]


df_signal = DataFrame(
    Esum = [bin1 for (bin1, bin2) in bin1_bin2],
    phi = [bin2 for (bin1, bin2) in bin1_bin2],
    eff = counts
)
CSV.write(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/figs/for_miro_ND_signal_2Dhist_sumE_vs_phi.csv"), df_signal)

# do the same for background histogram
binedges1_b = midpoints(binedges(h_background)[1]) |> collect
binedges2_b = midpoints(binedges(h_background)[2]) |> collect
bin1_bin2_b = [ (bin_mid_1, bin_mid_2) for bin_mid_1 in binedges1_b for bin_mid_2 in binedges2_b ]
counts_b = [lookup(h_background, bin1, bin2) for (bin1, bin2) in bin1_bin2_b]
df_background = DataFrame(
    Esum = [bin1 for (bin1, bin2) in bin1_bin2_b],
    phi = [bin2 for (bin1, bin2) in bin1_bin2_b],
    bkg_counts = counts_b
)
CSV.write(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/figs/for_miro_ND_background_2Dhist_sumE_vs_phi.csv"), df_background)



function make_nd_histogram(p, binning, vars)
    extracted_data = map(v -> getproperty.(p.data, Symbol(v)) |> collect, Symbol.(vars))
    data_tuple = Tuple(extracted_data)
    h = StatsBase.fit(Histogram, data_tuple, binning)

    return h
end

# create a tuple of vectors for each variable by iterating events in data


binning = ( 
    # 0:1:180,
    2400:10:3500,
    # 0:5:150,
    # 0:5:150,
    0:0.2:10,
    0:0.5:50,
)

total_number_of_bins = prod(length.(binning) .- 1)

h_signal_nd = make_nd_histogram(background[3], binning, ["sumE", "lPint", "lPext"]) 

