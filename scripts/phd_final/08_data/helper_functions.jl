
using Measurements, Statistics

mutable struct fit_exp
    phase::Int
    data_hist::Hist1D
    fitted_hists::Vector{Hist1D}
    exp_n::Vector{Measurement{Float64}}
    components
end

mutable struct fit_result
    phase
    data_hist
    fitted_hists
    params
    fit_n
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
    smeared_energy = fwhm == zero(eltype(fwhm)) ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm, e_range)
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


function compute_pull_values(data_hist, fit_hist; min_expected = 1e-9)
    bincount_data = bincounts(data_hist)
    bincount_fit = bincounts(fit_hist)

    binerrors_data = binerrors(data_hist)
    binerrors_fit = binerrors(fit_hist)

    pulls = Float64[]

    for i in eachindex(bincount_data)
        O = bincount_data[i]
        E = bincount_fit[i]

        if E <= min_expected
            continue
        end

        sigma2 = binerrors_data[i]^2 + binerrors_fit[i]^2
        if sigma2 <= 0
            continue
        end

        push!(pulls, (O - E) / sqrt(sigma2))
    end

    return pulls
end


function plot_pull_distribution(;
    data_hist,
    fit_hist,
    outputpath = nothing,
    title = "Pull distribution",
    xlabel = "(Data - Fit) / σ",
    ylabel = "Bins",
    bins = (-5:0.5:5),
    figure_size = (950, 700),
    fontsize = 28,
    color = "#005C8A",
    show_standard_normal = true,
    px_per_unit = 3,
)
    pulls = compute_pull_values(data_hist, fit_hist)

    f = Figure(size = figure_size, fontsize = fontsize)
    ax = Axis(f[1, 1], title = title, xlabel = xlabel, ylabel = ylabel)

    h_pull = Hist1D(pulls; binedges = bins)
    CairoMakie.hist!(
        ax,
        h_pull,
        color = (color, 0.4),
        strokewidth = 2,
        strokecolor = color,
    )

    if show_standard_normal && length(pulls) > 0
        x = range(first(bins), last(bins), length = 800)
        binw = bins[2] - bins[1]
        y = pdf.(Normal(0, 1), x) .* length(pulls) .* binw
        lines!(ax, x, y, color = :black, linewidth = 3, linestyle = :dash)
    end

    if length(pulls) > 0
        pull_mean = mean(pulls)
        pull_std = std(pulls)
        text!(
            ax,
            0.03,
            0.95,
            space = :relative,
            align = (:left, :top),
            text = "N = $(length(pulls))\nmean = $(round(pull_mean, digits = 2))\nσ = $(round(pull_std, digits = 2))",
            fontsize = fontsize - 4,
        )
    end

    CairoMakie.vlines!(ax, [0.0], color = :black, linestyle = :dot)

    if !isnothing(outputpath)
        save(outputpath, f, px_per_unit = px_per_unit)
    end

    return f, pulls
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
    data_markersize = 14,
    data_whiskerwidth = 8,
    component_linewidth = 4,
    total_fit_color = :red,
    total_fit_linewidth = 6,
    stacked = false,
    logy = false,
    filled = false,
    legend_patchsize = (30, 40, 100, 0),
    legend_rowgap = 20,
    legend_tellheight = false,
    legend_tellwidth= true,
    ratio_reference = 1.0,
    legend_width = 750,
    legend_nbanks = 1,
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

        el = filled ? PolyElement(color = c, strokecolor = c) : LineElement(color = c, linewidth = component_linewidth+5)
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
            "χ² / ndf \n= ",
            round(chi_ndf[2], digits = 2),
            " / ",
            chi_ndf[3],
            "\n= ",
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
        LineElement(color = total_fit_color, linewidth = total_fit_linewidth+5)

    legend_summary_entries = Any[
        MarkerElement(color = :black, marker = :circle, markersize = data_markersize+5),
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
    f[1:2, 2] = legend_grid

    Legend(legend_grid[1, 1], legend_summary_entries, legend_summary_labels,"Summary",
        tellheight = legend_tellheight,
        patchsize = legend_patchsize,
        rowgap = legend_rowgap,
        tellwidth = legend_tellwidth,
        patchlabelgap = 20,
        width = legend_width

    )

    Legend(legend_grid[2, 1], legend_entries, legend_labels,"Components",
        tellheight = legend_tellheight,
        patchsize = legend_patchsize,
        rowgap = legend_rowgap,
        tellwidth = legend_tellwidth,
        patchlabelgap = 20,
        width = legend_width,
        nbanks = legend_nbanks,
    )

    rowsize!(f.layout, 2, Relative(0.25))
    colgap!(f.layout, 1, Relative(0.05))

    save(outputpath, f, px_per_unit = px_per_unit)

    return f
end


function plot_fit_pull(;
    phase,
    data_hist,
    fitted_hists,
    outputpath,
    component_labels,
    title = nothing,
    xlabel = "Sum Energy (keV)",
    ylabel = "Counts / bin",
    residual_ylabel = "Data - Fit",
    pull_ylabel = "Pull",
    main_limits = (0, 4000, 0, nothing),
    residual_limits = nothing,
    pull_limits = (-5.0, 5.0),
    blinded_roi = (2700, 3000),
    figure_size = (1850, 1700),
    fontsize = 30,
    colors = ["#0F4C81", "#CC7A00", "#2F4858", "#3D8C40", "#7A2E8C", "#A23E2A", "#00798C"],
    total_fit_color = :black,
    total_fit_linewidth = 5,
    component_linewidth = 3,
    data_markersize = 12,
    data_whiskerwidth = 7,
    residual_band_color = (:gray70, 0.4),
    pull_band_color = (:gray80, 0.45),
    legend_position = :rt,
    stacked = false,
    filled = false,
    logy = false,
    show_chi2 = true,
    chi2_nparams = nothing,
    chi2_text_pos = (0.76, 0.90),
    show_pull_stats = true,
    px_per_unit = 3,
)

    use_stacked = stacked
    total_fit = sum(fitted_hists)

    x_edges = binedges(data_hist)
    x_step = x_edges[2] - x_edges[1]
    x_centers = collect(range(x_edges[1] + x_step / 2, x_edges[end] - x_step / 2, step = x_step))

    data_counts = bincounts(data_hist)
    fit_counts = bincounts(total_fit)
    data_err = binerrors(data_hist)
    fit_err = binerrors(total_fit)
    sigma_combined = sqrt.(data_err .^ 2 .+ fit_err .^ 2)

    residuals = data_counts .- fit_counts
    residual_err = sigma_combined
    pulls = [sigma_combined[i] > 0 ? residuals[i] / sigma_combined[i] : 0.0 for i in eachindex(residuals)]

    finite_residuals = filter(isfinite, residuals)
    residual_span = if isnothing(residual_limits)
        if isempty(finite_residuals)
            1.0
        else
            max(abs(minimum(finite_residuals)), abs(maximum(finite_residuals)), 1.0)
        end
    else
        nothing
    end

    resolved_residual_limits = isnothing(residual_limits) ?
        (x_edges[1], x_edges[end], -1.25 * residual_span, 1.25 * residual_span) :
        residual_limits

    main_axis_limits = if logy && !isnothing(main_limits[3]) && main_limits[3] <= 0
        (main_limits[1], main_limits[2], 0.1, main_limits[4])
    else
        main_limits
    end

    f = Figure(size = figure_size, fontsize = fontsize)

    ax_main = Axis(
        f[1, 1],
        ylabel = ylabel,
        title = isnothing(title) ? "Phase $phase" : title,
        limits = main_axis_limits,
        yscale = logy ? log10 : identity,
    )

    ax_res = Axis(
        f[2, 1],
        ylabel = residual_ylabel,
        limits = resolved_residual_limits,
    )

    ax_pull = Axis(
        f[3, 1],
        xlabel = xlabel,
        ylabel = pull_ylabel,
        limits = (x_edges[1], x_edges[end], pull_limits[1], pull_limits[2]),
    )

    plot_component_hist! = function (hist, color)
        if filled
            CairoMakie.hist!(
                ax_main,
                hist,
                color = (color, 0.35),
                strokewidth = component_linewidth,
                strokecolor = color,
            )
        else
            CairoMakie.stephist!(
                ax_main,
                hist,
                color = color,
                linewidth = component_linewidth,
            )
        end
    end

    legend_entries = Any[]
    legend_labels = Any[]

    for i in eachindex(fitted_hists)
        c = colors[mod1(i, length(colors))]
        hist_to_plot = use_stacked ? sum(fitted_hists[i:end]) : fitted_hists[i]
        plot_component_hist!(hist_to_plot, c)

        el = filled ? PolyElement(color = (c, 0.35), strokecolor = c) : LineElement(color = c, linewidth = component_linewidth + 2)
        push!(legend_entries, el)
        push!(legend_labels, component_labels[i])
    end

    if use_stacked
        CairoMakie.hist!(
            ax_main,
            total_fit,
            color = (total_fit_color, 0.16),
            strokewidth = 0,
        )
    else
        CairoMakie.stephist!(
            ax_main,
            total_fit,
            color = total_fit_color,
            linewidth = total_fit_linewidth,
        )
    end

    CairoMakie.scatter!(ax_main, data_hist, color = :black, markersize = data_markersize)
    CairoMakie.errorbars!(ax_main, data_hist, color = :black, whiskerwidth = data_whiskerwidth)

    CairoMakie.vspan!(ax_main, [blinded_roi[1]], [blinded_roi[2]], color = (:black, 0.15))
    CairoMakie.vspan!(ax_res, [blinded_roi[1]], [blinded_roi[2]], color = (:black, 0.15))
    CairoMakie.vspan!(ax_pull, [blinded_roi[1]], [blinded_roi[2]], color = (:black, 0.15))

    CairoMakie.band!(ax_res, x_centers, -residual_err, residual_err, color = residual_band_color)
    CairoMakie.scatter!(ax_res, x_centers, residuals, color = :black, markersize = data_markersize - 2)
    CairoMakie.errorbars!(ax_res, x_centers, residuals, residual_err, color = :black, whiskerwidth = data_whiskerwidth - 1)
    CairoMakie.hlines!(ax_res, [0.0], color = :black, linestyle = :dash, linewidth = 2)

    CairoMakie.band!(ax_pull, x_centers, fill(-1.0, length(x_centers)), fill(1.0, length(x_centers)), color = pull_band_color)
    CairoMakie.scatter!(ax_pull, x_centers, pulls, color = :black, markersize = data_markersize - 2)
    CairoMakie.errorbars!(ax_pull, x_centers, pulls, fill(1.0, length(x_centers)), color = :black, whiskerwidth = data_whiskerwidth - 1)
    CairoMakie.hlines!(ax_pull, [0.0], color = :black, linewidth = 2)
    CairoMakie.hlines!(ax_pull, [-2.0, 2.0], color = :gray35, linestyle = :dash, linewidth = 1.5)

    if show_chi2
        nparams = isnothing(chi2_nparams) ? length(fitted_hists) : chi2_nparams
        chi_ndf = get_chi2_ndf(data_hist, total_fit, nparams)
        chi_text = string(
            "χ² / ndf = ",
            round(chi_ndf[2], digits = 1),
            " / ",
            chi_ndf[3],
            " = ",
            round(chi_ndf[1], digits = 2),
        )
        text!(
            ax_main,
            chi2_text_pos[1],
            chi2_text_pos[2],
            text = chi_text,
            space = :relative,
            fontsize = fontsize - 4,
        )
    end

    if show_pull_stats
        pull_mean = mean(pulls)
        pull_std = std(pulls)
        text!(
            ax_pull,
            0.02,
            0.95,
            text = "μ = $(round(pull_mean, digits = 2)),  σ = $(round(pull_std, digits = 2))",
            space = :relative,
            align = (:left, :top),
            fontsize = fontsize - 6,
        )
    end

    summary_entries = Any[
        MarkerElement(color = :black, marker = :circle, markersize = data_markersize + 2),
        LineElement(color = total_fit_color, linewidth = total_fit_linewidth + 1),
        PolyElement(color = residual_band_color[1], strokecolor = :transparent),
        PolyElement(color = pull_band_color[1], strokecolor = :transparent),
    ]

    summary_labels = [
        "Data",
        "Total Fit",
        "Residual ±1σ",
        "Pull ±1σ",
    ]

    leg = Legend(
        f[1, 2],
        vcat(summary_entries, legend_entries),
        vcat(summary_labels, component_labels),
        "Model Components",
        patchsize = (32, 36, 80, 0),
        rowgap = 10,
        labelsize = fontsize - 8,
        titlesize = fontsize - 5,
        orientation = :vertical,
        tellheight = false,
        tellwidth = true,
        halign = :left,
        valign = :top,
    )
    # leg.position = legend_position

    linkxaxes!(ax_main, ax_res, ax_pull)
    hidexdecorations!(ax_main, grid = false)
    hidexdecorations!(ax_res, grid = false)

    rowgap!(f.layout, 10)
    rowsize!(f.layout, 1, Relative(0.58))
    rowsize!(f.layout, 2, Relative(0.22))
    rowsize!(f.layout, 3, Relative(0.20))
    colgap!(f.layout, 1, Relative(0.04))

    save(outputpath, f, px_per_unit = px_per_unit)
    return f
end






