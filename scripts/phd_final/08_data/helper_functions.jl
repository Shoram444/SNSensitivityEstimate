
using Measurements

mutable struct fit_result
    phase::Int
    data_hist::Hist1D
    fitted_hists::Vector{Hist1D}
    params::NamedTuple
    fit_n::Vector{Measurement{Float64}}
    components
end


function filter_data!(d, roi)
    d_new = deepcopy(d)
    for key in keys(roi)
        lo, hi = roi[key]
        filter!(row -> lo <= row[key] <= hi, d_new)
    end
    return d_new
end

function smear_energy(E::Real, fwhm::Real) 
    if(fwhm == 0)
        return E
    end
    sigma = fwhm/2.355 * E * sqrt(1000/E)
    return rand(Normal(E, sigma))
end

function get_smeared_energy(p, var, fwhm, e_range = (250, 5000))
    e1 = smear_energy.(getproperty.(p.data, :simulatedEnergy1), fwhm)
    e2 = smear_energy.(getproperty.(p.data, :simulatedEnergy2), fwhm)

    filter!(x -> e_range[1] <= x <= e_range[2], e1)
    filter!(x -> e_range[1] <= x <= e_range[2], e2)

    len1 = length(e1)
    len2 = length(e2)
    min_len = min(len1, len2)
    e1 = e1[1:min_len]
    e2 = e2[1:min_len]

    if var == :simulatedEnergy1 || var == :simulatedEnergy2
        return var == :simulatedEnergy1 ? e1 : e2
    elseif var == :sumEsimu
        return e1 .+ e2
    else
        return getproperty.(p.data, var)
    end
end

function get_expected_simu_counts(p, var, fwhm; bins = (0:100:4000), timeMeas = p.timeMeas, e_range = (250, 35000))
    smeared_energy = fwhm == 0 ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm, e_range)
    h = normalize(Hist1D(smeared_energy; binedges = bins), width = false)
    n_exp = p.activity * timeMeas * p.amount * (length(smeared_energy) / p.nTotalSim)
    return h * n_exp
end

function get_simu_histo(p, var, fwhm; bins = (0:100:4000), e_range = (250, 35000))
    smeared_energy = fwhm == 0 ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm, e_range)
    h = Hist1D(smeared_energy; binedges = bins)
    return h 
end

function get_chi2_ndf(data, fit, nparams)
    bincount_data = bincounts(data)
    bincount_fit = bincounts(fit)

    binerrors_data = binerrors(data)
    binerrors_fit = binerrors(fit)

    chi2 = 0.0
    ndf = 0

    for i in eachindex(bincount_data)
        O = bincount_data[i]
        E = bincount_fit[i]

        s = binerrors_data[i]^2 + binerrors_fit[i]^2

        if E > 5          # avoid division by zero
            chi2 += (O - E)^2 / s
            ndf += 1
        end
    end

    ndf = ndf - nparams
    return chi2 / ndf, chi2, ndf
end

function get_combined_bkg_hists(processes, binning)

    h_internal = Hist1D(;binedges = binning)
    h_Bi210 = Hist1D(;binedges = binning)
    h_detector = Hist1D(;binedges = binning)
    h_external = Hist1D(;binedges = binning)

    internals = ["Bi214_foil_bulk", "Tl208_foil_bulk"]
    Bi210 = ["Bi210_wire_surface"]
    detectors = ["Bi214_PMT_glass_bulk","Bi214_PMT5_glass_bulk","K40_PMT_glass_bulk","K40_PMT5_glass_bulk","Tl208_PMT_glass_bulk","Tl208_PMT5_glass_bulk","Bi214_calo_8inch_scin_bulk","K40_calo_8inch_scin_bulk", ]

    for p in processes
        h = Hist1D(p.dataVector; binedges = binning)

        if p.isotopeName in internals && integral(h) >0
            h_internal = h_internal + h
        elseif p.isotopeName in Bi210 && integral(h) >0
            h_Bi210 = h_Bi210 + h
        elseif p.isotopeName in detectors && integral(h) >0
            h_detector = h_detector + h
        elseif integral(h) >0
            h_external = h_external + h
        end
    end
    return vcat(h_internal, h_Bi210, h_detector, h_external)
end


function plot_lognormal_prior_comparison(;
    mean,
    relerr,
    outputpath = nothing,
    title = nothing,
    show_underlying_normal = false,
    lognormal_title = "Prior distribution",
    lognormal_xlabel = "x",
    ylabel = "Density",
    figure_size = show_underlying_normal ? (1800, 700) : (950, 700),
    fontsize = 28,
    color_lognormal = "#9A3A06",
    band_color_lognormal = (color_lognormal, 0.18),
    n_sigma = 4,
    n_points = 1200,
    px_per_unit = 3,
)

    σ_log = sqrt(log(1 + relerr^2))
    μ_log = log(mean) - 0.5 * σ_log^2

    normal_dist = Normal(μ_log, σ_log)
    lognormal_dist = LogNormal(μ_log, σ_log)


    lo_lognorm = max(eps(), quantile(lognormal_dist, 1e-4))
    hi_lognorm = quantile(lognormal_dist, 1 - 1e-4)
    x_lognorm = range(lo_lognorm, hi_lognorm, length = n_points)
    y_lognorm = pdf.(lognormal_dist, x_lognorm)

    lo68_lognorm, hi68_lognorm = quantile(lognormal_dist, [0.16, 0.84])
    lo95_lognorm, hi95_lognorm = quantile(lognormal_dist, [0.025, 0.975])

    fig = Figure(size = figure_size, fontsize = fontsize, figure_padding = 25)

    ax_lognorm = Axis(
        fig[1, 1],
        title = lognormal_title,
        xlabel = lognormal_xlabel,
        ylabel = ylabel,
    )

    if !isnothing(title)
        Label(fig[0, 1], title, fontsize = fontsize + 4, tellwidth = false)
    end

    band!(ax_lognorm, x_lognorm, zero.(y_lognorm), y_lognorm, color = band_color_lognormal)
    lines!(ax_lognorm, x_lognorm, y_lognorm, color = color_lognormal, linewidth = 4)
    vlines!(ax_lognorm, [mean], color = :black, linestyle = :dash, linewidth = 3)
    vlines!(ax_lognorm, [lo68_lognorm, hi68_lognorm], color = color_lognormal, linestyle = :dot, linewidth = 2)
    vlines!(ax_lognorm, [lo95_lognorm, hi95_lognorm], color = color_lognormal, linestyle = :dashdot, linewidth = 2)
    text!(
        ax_lognorm,
        0.03,
        0.95,
        space = :relative,
        align = (:left, :top),
        text = "mean = $(round(mean, sigdigits = 4))\nrel. err. = $(round(100 * relerr, digits = 1))%",
        fontsize = fontsize - 4,
    )

    axislegend(
        ax_lognorm,
        [
            LineElement(color = color_lognormal, linewidth = 4),
            LineElement(color = :black, linestyle = :dash, linewidth = 3),
            LineElement(color = color_lognormal, linestyle = :dot, linewidth = 2),
            LineElement(color = color_lognormal, linestyle = :dashdot, linewidth = 2),
        ],
        ["prior density", "input mean", "68% interval", "95% interval"],
        position = :rt,
    )

    if show_underlying_normal
        colgap!(fig.layout, 30)
    end

    if !isnothing(outputpath)
        save(outputpath, fig, px_per_unit = px_per_unit)
    end

    return fig, normal_dist, lognormal_dist
end




function plot_fit(;
    phase,
    data_hist,
    fitted_hists,
    fit_params,
    outputpath,
    component_labels,
    component_count_uncertainties = nothing,
    show_component_uncertainties = false,
    main_limits = (0, 4000, 0, nothing),
    ratio_limits = (0, 4000, 0.2, 1.8),
    blinded_roi = (2700, 3000),
    title = nothing,
    xlabel = "",
    ylabel = "Counts / 100 keV",
    ratio_xlabel = "Sum Energy (keV)",
    ratio_ylabel = "Data / Fit",
    figure_size = (2000, 1400),
    fontsize = 28,
    colors = ["#003865", "#FFB948", "#52473B", "#9A3A06", "#951272", "#006630", "#005C8A", "#FFB948"],
    px_per_unit = 3,
    chi2_nparams = nothing,
    show_chi2 = true,
    chi2_text_pos = (0.77, 0.9),
    include_component_counts = true,
    data_markersize = 10,
    data_whiskerwidth = 8,
    component_linewidth = 3,
    total_fit_color = :red,
    total_fit_linewidth = 4,
    stacked = false,
    logy = false,
    filled = false,
    legend_patchsize = (30, 40, 100, 0),
    legend_rowgap = 20,
    legend_tellheight = false,
    legend_tellwidth= true,
    ratio_reference = 1.0,
)

    use_stacked = stacked 
    total_fit = sum(fitted_hists)
    n_data = sum(bincounts(data_hist))
    n_mc_total = sum(bincounts(total_fit))
    main_axis_limits = if logy && !isnothing(main_limits[3]) && main_limits[3] <= 0
        (main_limits[1], main_limits[2], 0.1, main_limits[4])
    else
        main_limits
    end

    f = Figure(size = figure_size, fontsize = fontsize)

    ax = Axis(f[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = isnothing(title) ? "Phase $phase" : title,
        limits = main_axis_limits,
        yscale = logy ? log10 : identity,
    )

    ax_ratio = Axis(f[2, 1],
        xlabel = ratio_xlabel,
        ylabel = ratio_ylabel,
        limits = ratio_limits,
    )

    legend_entries = []
    legend_labels = []
    component_count_vals = nothing
    component_count_unc = nothing
    if show_component_uncertainties && !isnothing(component_count_uncertainties)
        raw_component_counts = collect(component_count_uncertainties)
        if all(x -> x isa Measurement, raw_component_counts)
            component_count_vals = Measurements.value.(raw_component_counts)
            component_count_unc = Measurements.uncertainty.(raw_component_counts)
        else
            component_count_unc = Float64.(raw_component_counts)
        end
    end

    plot_component_hist! = function (hist, color)
        if filled
            CairoMakie.hist!(
                ax,
                hist,
                color = color,
                strokewidth = component_linewidth,
                strokecolor = color,
            )
        else
            CairoMakie.stephist!(
                ax,
                hist,
                color = color,
                linewidth = component_linewidth,
            )
        end
    end

    for i in eachindex(fitted_hists)
        c = colors[mod1(i, length(colors))]

        hist_to_plot = use_stacked ? sum(fitted_hists[i:end]) : fitted_hists[i]
        plot_component_hist!(hist_to_plot, c)

        n_fit = isnothing(component_count_vals) ?
            sum(bincounts(fitted_hists[i])) :
            component_count_vals[i]
        lbl = if include_component_counts
            if !isnothing(component_count_unc) && i <= length(component_count_unc)
                rich(component_labels[i], "\n \n n = $(round(n_fit, digits = 1)) ± $(round(component_count_unc[i], digits = 1))")
            else
                rich(component_labels[i], "\n \n n = $(round(n_fit, digits = 1))")
            end
        else
            component_labels[i]
        end

        el = filled ? PolyElement(color = c, strokecolor = c) : LineElement(color = c, linewidth = component_linewidth)
        push!(legend_entries, el)
        push!(legend_labels, lbl)
    end

    if use_stacked
        CairoMakie.hist!(
            ax,
            total_fit,
            color = (total_fit_color, 0.15),
            strokewidth = 0,
            label = "Total Fit",
        )
    else
        CairoMakie.stephist!(
            ax,
            total_fit,
            color = total_fit_color,
            linewidth = total_fit_linewidth,
            label = "Total Fit",
        )
    end
    CairoMakie.scatter!(ax, data_hist,
        color = :black,
        markersize = data_markersize,
        label = "Data",
    )

    CairoMakie.errorbars!(ax, data_hist,
        color = :black,
        whiskerwidth = data_whiskerwidth,
    )
    CairoMakie.vspan!(ax, [blinded_roi[1]], [blinded_roi[2]], color = :black)
    CairoMakie.vspan!(ax_ratio, [blinded_roi[1]], [blinded_roi[2]], color = :black)

    if show_chi2
        nparams = isnothing(chi2_nparams) ?
            (fit_params isa NamedTuple ? length(keys(fit_params)) : length(fitted_hists)) :
            chi2_nparams

        chi_ndf = get_chi2_ndf(data_hist, total_fit, nparams)
        chi_text = string(
            "χ² / ndf = ",
            round(chi_ndf[2], digits = 2),
            " / ",
            chi_ndf[3],
            " = ",
            round(chi_ndf[1], digits = 2),
        )

        text!(ax, chi2_text_pos[1], chi2_text_pos[2],
            text = chi_text,
            space = :relative,
            fontsize = fontsize - 2,
        )
    end

    x_min = binedges(data_hist)[1]
    x_max = binedges(data_hist)[end]
    x_step = binedges(data_hist)[2] - binedges(data_hist)[1]

    xs = collect(range(x_min + x_step / 2, x_max - x_step / 2, step = x_step))
    ratio_hist = data_hist ./ total_fit
    ys = bincounts(ratio_hist)

    CairoMakie.scatter!(ax_ratio, xs, ys, color = :black)
    CairoMakie.errorbars!(
        ax_ratio,
        xs,
        ys,
        binerrors(ratio_hist),
        color = :black,
        whiskerwidth = data_whiskerwidth,
    )
    CairoMakie.hlines!(ax_ratio, [ratio_reference], linestyle = :dash, color = total_fit_color)

    total_legend_element = use_stacked ?
        PolyElement(color = (total_fit_color, 0.35), strokecolor = :transparent) :
        LineElement(color = total_fit_color, linewidth = total_fit_linewidth)

    legend_summary_entries = Any[
        MarkerElement(color = :black, marker = :circle, markersize = data_markersize),
        total_legend_element,
    ]
    n_mc_total_disp = isnothing(component_count_vals) ? n_mc_total : sum(component_count_vals)
    total_unc = isnothing(component_count_unc) ? nothing : sqrt(sum(component_count_unc .^ 2))
    legend_summary_labels = [
        "data, n = $(round(n_data, digits = 1))",
        isnothing(total_unc) ?
            "MC total, n = $(round(n_mc_total_disp, digits = 1))" :
            "MC total, n = $(round(n_mc_total_disp, digits = 1)) ± $(round(total_unc, digits = 1))",
    ]

    legend_grid = GridLayout()
    f[1, 2] = legend_grid

    Legend(legend_grid[1, 1], legend_summary_entries, legend_summary_labels,
        tellheight = legend_tellheight,
        patchsize = legend_patchsize,
        rowgap = legend_rowgap,
        tellwidth = legend_tellwidth,

    )

    Legend(legend_grid[2, 1], legend_entries, legend_labels,
        tellheight = legend_tellheight,
        patchsize = legend_patchsize,
        rowgap = legend_rowgap,
        tellwidth = legend_tellwidth,
        width = 350
    )

    rowsize!(f.layout, 2, Relative(0.25))
    save(outputpath, f, px_per_unit = px_per_unit)

    return f
end






