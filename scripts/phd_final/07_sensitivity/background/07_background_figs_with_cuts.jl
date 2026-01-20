###############################
# plot_background_spectra.jl #
###############################

using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate
using CairoMakie
using FHist
using ColorSchemes
using StatsBase

###############################
# Configuration helpers
###############################

STACK_COLORS = ColorSchemes.tableau_20

make_bins((low, high, bw)) = low:bw:(high + bw)


topology_roi = Dict(
    :dz => (0, 10000),
    :dy => (0, 10000),
    :sumE => (0, 10000),
    :phi => (0, 180),
    :Pint => (0, 1),
    :Pext => (0, 1),
    :trackLength1 => (0, 10000),
    :trackLength2 => (0, 10000),
)

tof_roi = Dict(
    :r => (0,100000), 
    :sumE => (0, 10000),
    :phi => (0, 180),
    :Pint => (0.04, 1),
    :Pext => (0, 0.01),
    :trackLength1 => (0, 150000),
    :trackLength2 => (0, 150000),
)

vertex_tof_roi = Dict(
    :r => (0,100), 
    :sumE => (0, 10000),
    :phi => (0, 180),
    :Pint => (0.04, 1),
    :Pext => (0, 0.01),
    :trackLength1 => (0, 150000),
    :trackLength2 => (0, 150000),
)

vertex_energy_tof_roi = Dict(
    :r => (0,100), 
    :sumE => (0, 3500),
    :phi => (0, 180),
    :Pint => (0.04, 1),
    :Pext => (0, 0.01),
    :trackLength1 => (0, 150000),
    :trackLength2 => (0, 150000),
)

vertex_lenfth_energy_tof_roi = Dict(
    :r => (0,100), 
    :sumE => (0, 3500),
    :phi => (0, 180),
    :Pint => (0.04, 1),
    :Pext => (0, 0.01),
    :trackLength1 => (0, 1500),
    :trackLength2 => (0, 1500),
)


nemo3_roi = Dict(
    :r => (0,50), 
    :sumE => (0, 3500),
    :phi => (0, 180),
    :Pint => (0.04, 1),
    :Pext => (0, 0.01),
    :trackLength1 => (0, 1500),
    :trackLength2 => (0, 1500),
)


###############################
# Data loading
###############################

function load_background_hists(
    simdir;
    mode::Symbol,
    roi::Dict,
    backgrounds::Vector{String},
    binning = (0, 4000, 100),
)

    low, high, bw = binning
    bins = make_bins(binning)

    all_processes = load_data_processes(
        simdir,
        String(mode);
        fwhm = 0.0,
        roi  = roi,
    )

    bkgs = [get_process(b, all_processes) |> first for b in backgrounds]

    for b in bkgs
        set_bins!(b, bins)
    end

    hists = get_bkg_counts_1D.(bkgs)
    for (i,h) in enumerate(hists)
        if isnan(integral(h))
            replace!(hists[i].bincounts, NaN => 0.0)
        end
    end


    min_cts = find_min_bincounts(hists)
    max_cts = 2.2 * find_max_bincounts(hists)

    return hists, min_cts, max_cts, bw
end

###############################
# Plotting helpers
###############################

function find_min_bincounts(hists)
    min_ = 1e200
    for hist in hists
        for count in bincounts(hist)
            if !isnan(count) && count > 0
                if min_ == 0 || count < min_
                    min_ = count
                end
            end
        end
    end
    return min_
end

function find_max_bincounts(hists::Hist1D)
    if typeof(hists) <: FHist.Hist1D
        return maximum(bincounts(hists))
    end
end

function find_max_bincounts(hists::Vector{Hist1D})
    max_ = 0
    for hist in hists
        for count in bincounts(hist)
            if !isnan(count) && count >= 0
                if max_ == 0 || count > max_
                    max_ = count
                end
            end
        end
    end
    return max_
end

function stacked_hist!(ax, hists; labels, colors, fill_area)
    if fill_area
        hist!(ax, sum(hists),
            label = isnothing(labels) ? nothing : labels[1],
            color = colors[1],
            strokewidth = 2,
            strokecolor = :black)
    else
        stephist!(ax, sum(hists),
            label = isnothing(labels) ? nothing : labels[1],
            color = colors[1],
            linewidth = 2)        
    end
    # hist!(ax, sum(hists), color = (:black, 0.10), )

    errorbars!(ax, sum(hists); whiskerwidth=7, color = colors[1], clamp_bincounts=true)
    # edges = binedges(hists[1]) |> midpoints
    # errs = binerrors(hists[1])
    # replace!(errs, NaN => 0.0)
    # replace!(errs, Inf => 0.0)
    # bcs = FHist.bincounts(sum(hists))
    # errorbars!(ax, edges, bcs, errs; whiskerwidth=7, color = colors[1])

    for i in 2:length(hists)
        if fill_area
            hist!(ax, sum(hists[i:end]),
                label = isnothing(labels) ? nothing : labels[i],
                color = colors[i],
                strokewidth = 2,
                strokecolor = :black)
        else
            stephist!(ax, sum(hists[i:end]),
                label = isnothing(labels) ? nothing : labels[i],
                color = colors[i],
                linewidth = 2)
        end
        nonzero_hists = sum(filter(h -> maximum(bincounts(h)) > 0, hists[i:end])) # when histogram is empty gives NaN errors
        errorbars!(ax, nonzero_hists; whiskerwidth=7, color = colors[i], clamp_bincounts=true)
        # errorbars!(ax, sum(hists[i:end]); whiskerwidth=7, color = colors[i], clamp_errors=false)
        # edges = binedges(hists[i]) |> midpoints
        # errs = abs.(binerrors(hists[i]))
        # replace!(errs, NaN => 0.0)
        # replace!(errs, Inf => 0.0)
        # @show i, errs

        # bcs = FHist.bincounts(sum(hists[i:end]))
        # errorbars!(ax, edges, bcs, errs; whiskerwidth=7, color = colors[i])

    end
end

function all_individual_hist!(ax, hists; labels, colors, fill_area)
    for i in 1:length(hists)
        if fill_area
            hist!(ax, hists[i],
                label = isnothing(labels) ? nothing : labels[i],
                color = colors[i],
                strokewidth = 1,
                strokecolor = :black)
        else
            stephist!(ax, hists[i],
                label = isnothing(labels) ? nothing : labels[i],
                color = colors[i],
                linewidth = 4)
        end
        errorbars!(ax, hists[i]; whiskerwidth=7, color = colors[i])

    end
end


function single_hist!(ax, hist; label, color)
    hist!(ax, hist,
          label = label,
          color = color,
          strokewidth = 1,
          strokecolor = :black)
end

function axis_labels(mode::Symbol)
    if mode == :phi
        return L"$\varphi ~(^{\circ})$", L"Counts $$"
    elseif mode == :singleE
        return L"Single electron energy (keV)$$", L"Counts $$"
    else
        return L"Summed 2-electron energy (keV)$$", L"Counts $$"
    end
end

###############################
# Plot generators
###############################

function plot_stacked(
    hists;
    labels,
    bw,
    min_cts,
    max_cts,
    mode::Symbol,
    scale::Symbol,
    colors,
    fill_area
)

    xlabel, ylabel = axis_labels(mode)

    f = Figure(size = (1400, 900), fontsize = 24, figure_padding = 24)
    ax = Axis(
        f[1, 1];
        xlabel,
        ylabel = L"%$ylabel/ (\mathrm{yr} \cdot %$bw \mathrm{keV})",
        title  = "Simulated SuperNEMO background",
        yscale = scale == :log ? log10 : identity,
    )

    stacked_hist!(ax, hists; labels, colors = colors, fill_area = fill_area)

    ylims!(ax, min_cts, scale == :log ? 1e4 * max_cts : max_cts)

    axislegend(ax, position = :lt, orientation = :horizontal, nbanks = 5)


    return f
end

function plot_individual(
    hist;
    label,
    bw,
    mode::Symbol,
    scale::Symbol,
    color,
)

    xlabel, ylabel = axis_labels(mode)

    f = Figure(size = (1200, 800), fontsize = 22, figure_padding = 24)
    ax = Axis(
        f[1, 1];
        xlabel,
        ylabel = L"%$ylabel/ (\mathrm{yr} \cdot %$bw \mathrm{keV})",
        title  = label,
        yscale = scale == :log ? log10 : identity,
    )

    single_hist!(ax, hist; label, color)

    min_cts, max_cts = find_min_bincounts([hist]), find_max_bincounts([hist])

    if scale == :log
        ylims!(ax, max(min_cts, 1e-1), 1e2 * max_cts)
    end

    axislegend(ax, position = :rt)

    return f
end

###############################
# Main user-facing function
###############################

function make_background_spectra(
    simdir;
    roi::Dict,
    analysis_dict::Dict,
    backgrounds::Vector{String},
    binning     = (0, 4000, 100),
    plot_mode   = :stacked,     # :stacked or :individual
    outdir      = pwd(),
    scale = :identity,
    fig_format = :png
)

    mode = Symbol(analysis_dict[:mode])
    plot_mode = Symbol(plot_mode)

    hists, min_cts, max_cts, bw =
        load_background_hists(
            simdir;
            mode,
            roi,
            backgrounds,
            binning,
        )

    save_tag = savename(analysis_dict)



    if plot_mode == :stacked
        fig = plot_stacked(
            hists;
            labels   = backgrounds,
            bw       = bw,
            min_cts  = min_cts,
            max_cts  = max_cts,
            mode     = mode,
            scale    = scale,
            colors   = STACK_COLORS,
        )
        
        if(fig_format == :svg)
            safesave(joinpath(outdir, "stacked","background_$(save_tag)_$(scale).svg"), fig)
        else
            safesave(joinpath(outdir, "stacked","background_$(save_tag)_$(scale).png"), fig)
        end
        
    elseif plot_mode == :individual
        for (i, hist) in enumerate(hists)
            fig = plot_individual(
                hist;
                label   = backgrounds[i],
                bw      = bw,
                # min_cts = min_cts,
                # max_cts = max_cts,
                mode    = mode,
                scale   = scale,
                color   = STACK_COLORS[i],
            )

            if(fig_format == :svg)
                safesave(joinpath(outdir, "individual",
                    "background_$(backgrounds[i])_$(save_tag)_$(scale).svg"), fig)
            else
                safesave(joinpath(outdir, "individual",
                    "background_$(backgrounds[i])_$(save_tag)_$(scale).png"), fig)
            end
        end
    else
        error("Unknown plot_mode = $plot_mode (use :stacked or :individual)")
    end
end

###############################
# Run 
###############################


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

    # external gamma
    "Tl208_hall_bulk", 
    "Bi214_hall_bulk",
    "K40_hall_bulk",
]



analysisDict    = Dict(:mode => "sumE", :cuts => "topology_roi")

if( analysisDict[:cuts] == "nemo3_roi" )
    roi = nemo3_roi
elseif( analysisDict[:cuts] == "topology_roi" )
    roi = topology_roi
elseif( analysisDict[:cuts] == "tof_roi" )
    roi = tof_roi
elseif( analysisDict[:cuts] == "vertex_tof_roi" )
    roi = vertex_tof_roi
elseif( analysisDict[:cuts] == "vertex_energy_tof_roi" )
    roi = vertex_energy_tof_roi
elseif( analysisDict[:cuts] == "vertex_length_energy_tof_roi" )
    roi = vertex_lenfth_energy_tof_roi
end


# make_background_spectra(
#     datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41");
#     roi         = roi,
#     analysis_dict = analysisDict,
#     backgrounds = backgrounds,
#     binning     = (0, 4000, 100),
#     plot_mode   = :stacked,     # :stacked ,:individual
#     outdir      = scriptsdir("phd_final/07_sensitivity/background/figs"),
#     scale       = :identity,
#     fig_format  = :png
# )



### final figure with categories
labels_latex = [
    # source foil
    L"2\nu\beta\beta", 
    L"^{214}\mathrm{Bi}", 
    L"^{40}\mathrm{K}", 
    L"^{234m}\mathrm{Pa}", 
    L"^{208}\mathrm{Tl}", 

    # radon
    L"^{210}\mathrm{Bi}",     
    L"^{214}\mathrm{Bi}", 

    # detector contamination
    L"\mathrm{PMT}",
    L"\mathrm{Scintillator}",

    # external gamma
    L"^{208}\mathrm{Tl}",
    L"^{214}\mathrm{Bi}", 
    L"^{40}\mathrm{K}",
]

categories = [
    "source foil",
    "radon",
    "detector contamination",
    "external gamma"
]

# let 
begin
    
    begin
        binning = (300, 3500, 100)
        # binning = (0, 180, 5)
        simdir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41")
        mode = analysisDict[:mode] 
        fill_area = false
        stacked = false
        # roi = roi
        scale = :log10
        cat_colors = ["#4E79A7","#C7E6E3","#9AD0CB","#76B7B2","#4E8F8B","#F28E2B","#C66A00","#E15759","#B22222", "#7BC36D", "#59A14F", "#3F7F3A"]
        # cat_colors = [STACK_COLORS[1], STACK_COLORS[2], STACK_COLORS[3]]
    end
    
    low, high, bw = binning
    # low, high, bw = 2800:100:3200
    bins = make_bins(binning) # 0:100:3200

    all_processes = load_data_processes(
        simdir,
        String(mode);
        fwhm = 0.0,
        roi  = roi,
    )

    bkgs = [get_process(b, all_processes) |> first for b in backgrounds]

    for b in bkgs
        set_bins!(b, bins)
    end

    tot = 0

    
    hists = get_bkg_counts_1D.(bkgs) .* inv(2.8)  # normalize to counts / yr
    for (i,h) in enumerate(hists)
        if isnan(integral(h))
            replace!(hists[i].bincounts, NaN => 0.0)
        end
        if occursin("tof", analysisDict[:cuts]) && i in (16:18)  # external gamma with tof cuts
            hists[i].bincounts .= hists[i].bincounts * inv(20)  # aritfitial scaling to ~500 ToF (this data gives only ~20 cause of low statistics)
        end 
        npass = length(bkgs[i].dataVector)
        np_err = measurement(sqrt(npass), sqrt(npass))
        nsim = bkgs[i].nTotalSim
        nsim_err = measurement(sqrt(nsim), sqrt(nsim))
        eff_err = npass/nsim

        
        measure_cts = measurement.(bincounts(h), binerrors(h))
        println(backgrounds[i], ", eff = ", eff_err, ", counts = ", sum(measure_cts))
        tot += sum(measure_cts)
    end
    println("Total cts = ", tot)


    for (i,h) in enumerate(hists)
        if i in (1:7)
            continue
        end
        if (roi != "topology_roi" && mode == "sumE")
            # zero out bins outside 2.8-3.2 MeV for better visibility
            for (j,edge) in enumerate(midpoints(binedges(h)))
                if edge > 2700
                    h.bincounts[j] = 0.0
                    h.sumw2[j] = 0.0
                end
            end
        end
    end

    sum_hists = vcat(
        hists[1:5],   # source foil
        hists[6:7],   # radon
        sum(hists[8:13]),  # detector contamination PMT
        sum(hists[14:15]),  # detector contamination Scintillator
        hists[16:18],  # external gamma
    )

    min_cts = max(1e-4,find_min_bincounts(sum_hists))
    max_cts = find_max_bincounts(sum(sum_hists))

    f = Figure(size = (1400, 900), fontsize = 36, figure_padding = 16)

    legend_foil = [
        [
            LineElement(color = cat_colors[i], linewidth = 3), 
            LineElement(color = cat_colors[i], points = Point2f[(0.5, 0), (0.5, 1)],linewidth = 3) 
        ] for i in 1:5
    ] 


    # legend_foil = [PolyElement(polycolor = cat_colors[i]) for i in 1:5] 
    legend_radon = [
        [
            LineElement(color = cat_colors[i], linewidth = 3), 
            LineElement(color = cat_colors[i], points = Point2f[(0.5, 0), (0.5, 1)],linewidth = 3) 
        ] for i in 6:7
    ]
    legend_detector = [
        [
            LineElement(color = cat_colors[i], linewidth = 3), 
            LineElement(color = cat_colors[i], points = Point2f[(0.5, 0), (0.5, 1)],linewidth = 3) 
        ] for i in 8:9
    ]
    legend_gamma = [
        [
            LineElement(color = cat_colors[i], linewidth = 3), 
            LineElement(color = cat_colors[i], points = Point2f[(0.5, 0), (0.5, 1)],linewidth = 3) 
        ] for i in 10:12
    ]



    legend_master = Legend(
        f,
        [
            legend_foil,
            legend_radon,
            legend_detector,
            legend_gamma,
        ],
        [
            labels_latex[1:5],
            labels_latex[6:7],
            labels_latex[8:9],
            labels_latex[10:12],
        ],
        [
            "Source foil",
            "Radon",
            "Detector contamination",
            "External gamma",
        ],
        patchsize = (45, 35),
        rowgap = 20,
        tellheight = true,
        padding = (10, 10, 10, 10)
    )

    legend_master.nbanks = 3
    legend_master.tellheight = false
    legend_master.tellwidth = true

    title = if analysisDict[:cuts] == "nemo3_roi"
        "SuperNEMO: Simulated background for 1 yr of exposure \nNEMO-3 like cuts"
    elseif analysisDict[:cuts] == "topology_roi"
        "SuperNEMO: Simulated background for 1 yr of exposure \ntopology only cuts"
    elseif analysisDict[:cuts] == "tof_roi"
        "SuperNEMO: Simulated background for 1 yr of exposure \ntopology + ToF cuts"
    elseif analysisDict[:cuts] == "vertex_tof_roi"
        "SuperNEMO: Simulated background for 1 yr of exposure \nvertex + ToF cuts"
    elseif analysisDict[:cuts] == "vertex_energy_tof_roi"
        "SuperNEMO: Simulated background for 1 yr of exposure \ntopology + ToF + vertex + energy cuts"
    elseif analysisDict[:cuts] == "vertex_length_energy_tof_roi"
        "SuperNEMO: Simulated background for 1 yr of exposure \nvertex + length + energy + ToF cuts"
    end

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
    axislegend(ax, position = :rt, padding = (20,20,20,20))

    if stacked
        stacked_hist!(ax, sum_hists; labels = nothing, colors = cat_colors, fill_area = fill_area)
    else
        all_individual_hist!(ax, sum_hists; labels = nothing, colors = cat_colors, fill_area = fill_area)
    end


    if scale == :log10
        ylims!(ax, 10^-4, 1e2 * max_cts)
        decades = -4:ceil(log10(1e2*max_cts)) .|> Int
        ticks   = 10.0 .^ decades 
        labels  = [L"10^{%$d}" for d in decades]
        ax.yticks = ( ticks,labels )
        

    else
        ylims!(ax, 0, 1.2*max_cts)  
    end

    xlims!(ax, binning[1], binning[2])
    # ax.xticks = binning[1]:binning[3]:binning[2]
    # ax2 = Axis(f, bbox=BBox(1400*0.75, 1400*0.95, 900*0.5,900*0.85), limits = (2600, 4000, 0, 10), title = "Zoomed-in view")
    # stacked_hist!(ax2, sum_hists; labels = nothing, colors = cat_colors)

    f[1, 2] = legend_master

    save(joinpath(scriptsdir("phd_final/07_sensitivity/background/figs/stacked"),
        "mode_$(mode)_background_$(analysisDict[:cuts])_categories_$(scale)_fill_$(fill_area)_stacked_$(stacked).png"), f)

    f
end


bkgs = [get_process(b, all_processes) |> first for b in backgrounds]

for b in bkgs
    set_bins!(b, bins)
end

hists = get_bkg_counts_1D.(bkgs) .* inv(2.8)  # normalize to counts / yr
for (i,h) in enumerate(hists)
    if isnan(integral(h))
        replace!(hists[i].bincounts, NaN => 0.0)
    end
end

sum_hists = vcat(
    hists[1:5],   # source foil
    hists[6:7],   # radon
    sum(hists[8:13]),  # detector contamination PMT
    sum(hists[14:15]),  # detector contamination Scintillator
)



f = Figure()
ax1 = Axis(f[1,1], yscale = log10, )
stacked_hist!(ax1, sum_hists; labels = nothing, colors = cat_colors, fill_area = fill_area)
# ylims!(ax1, min_cts, 1e2 * max_cts)
ylims!(ax1, 10^-2, 10^2)
decades = -2:2
ticks   = 10.0 .^ decades
labels  = [L"10^{%$d}" for d in decades]
ax1.yticks = ( ticks,labels )
f


ylims!(ax1, 10^-3, 1e2 * max_cts)
decades = -2:ceil(log10(1e2*max_cts)) .|> Int
ticks   = 10.0 .^ decades 
labels  = [L"10^{%$d}" for d in decades]
ax1.yticks = ( ticks,labels )

f


### 
nonzero_hist = sum_hists[1]
f = Figure()
ax = Axis(f[1,1], yscale = identity, )
# stephist!(ax, nonzero_hist)
errorbars!(ax, nonzero_hist; whiskerwidth=7, color = :red, clamp_bincounts = true, clamp_errors=true)

edges = binedges(nonzero_hist) |> midpoints
errs = binerrors(nonzero_hist)
bc = bincounts(sum(sum_hists[9:end]))
errorbars!(ax, edges, bc, errs, errs; whiskerwidth=7, color = :black)
stephist!(ax, sum(sum_hists[9:end]), color = :black, linewidth=2)
f

errorbars!(ax, nonzero_hist; whiskerwidth=7, color = :red, clamp_bincounts=true)
f


# errors on the total for cut set 1
m_set1 = measurement.(
    [
        5.54e4 ± 13.19,
        5.42±0.01,
        0.72±0.02,
        9.09e3±34.28,
        1.05e4±19.37,
        142.2±0.3,
        14.7±0.1,
        (6.82±0.07)*1e3,
        (3.01±0.04)*1e4,
        (1.52±0.02)*1e4,
        (1.37±0.01)*1e5,
        (3.98±0.05)*1e4,
        (4.62±0.06)*1e4,
        (2.35±0.05)*1e3,
        (4.99±0.19)*1e4,
        1.1±0.01,
        7.52±0.08,
        1.04±0.01

    ]
) |> sum

