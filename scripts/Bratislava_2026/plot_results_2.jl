using CSV
using DataFrames
using SNSensitivityEstimate
using CairoMakie
using LaTeXStrings
using Printf

base_dir = @__DIR__
sens_df = CSV.read(joinpath(base_dir, "nd_sensitivity.csv"), DataFrame)
world_bests_df = CSV.read(joinpath(base_dir, "world_bests.csv"), DataFrame)

times = collect(range(0.01, 5.0, 300))
t_base = SNparams["tYear"]
foil_mass = SNparams["foilMass"]
exposures = times .* foil_mass

# First 15 months of running at high-radon level.
switch_time_years = 15.0 / 12.0

# Radon levels in Bq/m^3 from nd_sensitivity.csv
const RADON_LOW = 0.00015  # 0.15 mBq/m^3
const RADON_MID = 0.00108  # 1.08 mBq/m^3
const RADON_HIGH = 0.03821 # 38.21 mBq/m^3

colors = ["#041E42", "#BE4D00", "#951272", "#006630", "#005C8A", "#FFB948", "#605643", "#302D23"]

function signal_tex(signal::AbstractString)
    if signal == "bb0nu_foil_bulk"
        return L"0\nu\beta\beta"
    elseif signal == "bb0nuM1_foil_bulk"
        return L"0\nu\beta\beta\chi^0"
    elseif signal == "bb0nuM2_foil_bulk"
        return L"0\nu\beta\beta\chi^0\chi^0"
    elseif signal == "RH050_foil_bulk" || signal == "RH040_foil_bulk"
        return L"\nu_L\nu_R\beta\beta"
    end
    return L"signal"
end

function title_label_tex(signal::AbstractString)
    if signal == "bb0nu_foil_bulk"
        return L"Simulated\ SuperNEMO\ sensitivity\ to\ \textrm{0\nu\beta\beta}\ decay $$"
    elseif signal == "bb0nuM1_foil_bulk"
        return L"Simulated\ SuperNEMO\ sensitivity\ to\ \textrm{0\nu\beta\beta\chi^0}\ decay $$"
    elseif signal == "bb0nuM2_foil_bulk"
        return L"Simulated\ SuperNEMO\ sensitivity\ to\ \textrm{0\nu\beta\beta\chi^0\chi^0}\ decay $$"
    elseif signal == "RH050_foil_bulk" || signal == "RH040_foil_bulk"
        return L"Simulated\ SuperNEMO\ sensitivity\ to\ \textrm{\nu_L\nu_R\beta\beta}\ decay $$"
    end
    return L"Simulated\ SuperNEMO\ sensitivity $$"
end

function world_best_label_tex(project::AbstractString, signal::AbstractString)
    if project == "CUPID" && signal == "bb0nu_foil_bulk"
        return L"$$\ CUPID-0 @ 90\%\ CI\ bayesian\ limit"
    elseif project == "CUPID" && signal == "bb0nuM1_foil_bulk"
        return L"$$\ CUPID\ world\ best\ @\ 90\%\ CI"
    elseif project == "CUPID" && signal == "bb0nuM2_foil_bulk"
        return L"$$\ CUPID\ world\ best\ @\ 90\%\ CI"
    end
    return L"$$\ world\ best\ @\ 90\%\ CI"
end

function legend_title_tex(signal::AbstractString)
    if signal == "bb0nu_foil_bulk"
        return L"SuperNEMO\ $$0\nu\beta\beta$$\ sensitivity\\ frequentist\ limit\ @\ 90\%\ CL"
    elseif signal == "bb0nuM1_foil_bulk"
        return L"SuperNEMO\ $$0\nu\beta\beta\chi^0$$\ sensitivity\\ frequentist\ limit\ @\ 90\%\ CL"
    elseif signal == "bb0nuM2_foil_bulk"
        return L"SuperNEMO\ $$0\nu\beta\beta\chi^0\chi^0$$\ sensitivity\\ frequentist\ limit\ @\ 90\%\ CL"
    elseif signal == "RH050_foil_bulk" || signal == "RH040_foil_bulk"
        return L"SuperNEMO\ $$\nu_L\nu_R\beta\beta$$\ sensitivity\\ frequentist\ limit\ @\ 90\%\ CL"
    end
    return L"SuperNEMO\ sensitivity"
end

function signal_filename(signal::AbstractString)
    return replace(signal, r"[^A-Za-z0-9]+" => "_")
end

function decade_ticks(y_low, y_high)
    p_low = Int(floor(log10(y_low)))
    p_high = Int(ceil(log10(y_high)))

    vals = [10.0^p for p in p_low:p_high]
    labels = [latexstring("10^{", string(p), "}") for p in p_low:p_high]
    return vals, labels
end

function y_limits_log(curves::Vector{Vector{Float64}})
    y_min = minimum(minimum(c) for c in curves)
    y_max = maximum(maximum(c) for c in curves)

    # Tight data-driven bounds.
    low_raw = y_min * 0.65
    high_raw = y_max * 1.12

    # Snap bounds to decade values and enforce at least 2 decades (factor 100).
    low = 10.0^floor(log10(low_raw))
    high = 10.0^ceil(log10(high_raw))
    if high / low < 100.0
        high = low * 100.0
    end

    return low, high
end

function get_row_for_radon(signal_rows::DataFrame, target_radon::Float64)
    idx = argmin(abs.(signal_rows.radon_level .- target_radon))
    return signal_rows[idx, :]
end

function cumulative_piecewise_background(
    bkg_first::Float64,
    bkg_after::Float64,
    t::Float64,
    switch_t::Float64,
    t_base::Float64,
)
    if t <= switch_t
        return bkg_first / t_base * t
    end
    b_first = bkg_first / t_base * switch_t
    b_after = bkg_after / t_base * (t - switch_t)
    return b_first + b_after
end

function compute_sensitivity_curve_piecewise(
    eff::Float64,
    bkg_first::Float64,
    bkg_after::Float64,
    times,
    t_base;
    switch_t = switch_time_years,
    efffactor = 1.0,
    extra_bkg = 0.0,
)
    α = 1.64485362695147
    t_halfs = Float64[]
    sizehint!(t_halfs, length(times))

    for t in times
        cumulative_bkg = cumulative_piecewise_background(bkg_first, bkg_after, t, switch_t, t_base)
        t_half = get_tHalf(
            SNparams["W"],
            SNparams["foilMass"],
            SNparams["Nₐ"],
            t,
            SNparams["a"],
            eff * efffactor,
            cumulative_bkg + extra_bkg,
            α;
            approximate = "table",
        )
        push!(t_halfs, t_half)
    end

    return t_halfs
end

function make_signal_plot(signal_rows::DataFrame, world_best_rows::DataFrame, signal::AbstractString)
    low = get_row_for_radon(signal_rows, RADON_LOW)
    mid = get_row_for_radon(signal_rows, RADON_MID)
    high = get_row_for_radon(signal_rows, RADON_HIGH)

    # Requested scenarios:
    # 1) 0.15 mBq only
    curve_low_only = compute_sensitivity_curve_piecewise(low.eff, low.bkgCounts, low.bkgCounts, times, t_base)

    # 2) Current: 38 mBq for first 15 months, then 1.08 mBq for the rest
    curve_current = compute_sensitivity_curve_piecewise(mid.eff, high.bkgCounts, mid.bkgCounts, times, t_base)

    # 3) Transition ideal: 38 mBq for first 15 months, then 0.15 mBq for the rest
    curve_ideal_transition = compute_sensitivity_curve_piecewise(low.eff, high.bkgCounts, low.bkgCounts, times, t_base)

    curves = [curve_low_only, curve_current, curve_ideal_transition]
    y_low, y_high = y_limits_log(curves)
    yt_vals, yt_labels = decade_ticks(y_low, y_high)

    fig = Figure(
        size = (950 * 1.3, 600),
        fontsize = 22,
        figure_padding = 20,
        fonts = (; regular = "TeX"),
    )

    ax = Axis(
        fig[1, 1],
        xlabel = L"$$\ exposure\ (\text{kg}\cdot\text{yr})",
        ylabel = L"Simulated\ sensitivity\ (yr) $$",
        yscale = log10,
        limits = (0, maximum(exposures), y_low, y_high),
        yminorticks = IntervalsBetween(10),
        yminorticksvisible = true,
        yminorticksize = 4,
    )

    ax_top = Axis(
        fig[1, 1],
        xlabel = L"$$time\ (yr)",
        xaxisposition = :top,
        yscale = log10,
        limits = (0, maximum(times), y_low, y_high),
    )

    Label(fig[0, 1], text = title_label_tex(signal))

    vline_positions = Float64[]
    if minimum(exposures) <= 17.5 <= maximum(exposures)
        push!(vline_positions, 17.5)
    end
    if minimum(exposures) <= 6.25 <= maximum(exposures)
        push!(vline_positions, 6.25)
    end
    if !isempty(vline_positions)
        vlines!(ax, vline_positions, linestyle = (:dash, :dense), color = (:black, 0.7), linewidth = 1.5)
        text!(17.5, y_low * 1.07, text = L"17.5\ \mathrm{kg}\cdot\mathrm{yr}", align = (:center, :top), rotation = pi / 2, fontsize = 18, color = (:black, 0.3))
        text!(6.25, y_low * 1.07, text = L"6.25\ \mathrm{kg}\cdot\mathrm{yr}", align = (:center, :top), rotation = pi / 2, fontsize = 18, color = (:black, 0.3))
    end

    p1 = lines!(ax, exposures, curve_low_only, color = colors[1], linewidth = 5.8)
    lines!(ax_top, times, curve_low_only, color = colors[1], linewidth = 5.8)

    p2 = lines!(ax, exposures, curve_current, color = colors[2], linewidth = 5.8, linestyle = (:dash, :dense))
    lines!(ax_top, times, curve_current, color = colors[2], linewidth = 5.8, linestyle = (:dash, :dense))

    p3 = lines!(ax, exposures, curve_ideal_transition, color = colors[3], linewidth = 5.8, linestyle = (:dashdot, :dense))
    lines!(ax_top, times, curve_ideal_transition, color = colors[3], linewidth = 5.8, linestyle = (:dashdot, :dense))

    handles = Any[p1, p2, p3]
    labels = Any[
        L"$$\ Rn\ @\ 0.15\ \mathrm{mBq/m^3}\ (all\ time)",
        L"$$\ 38.21\ \mathrm{mBq/m^3}\ for\ 15\ months\\ +\ 1.08\ \mathrm{mBq/m^3}\ after",
        L"$$\ 38.21\ \mathrm{mBq/m^3}\ for\ 15\ months\\ +\ 0.15\ \mathrm{mBq/m^3}\ after",
    ]

    for wb in eachrow(world_best_rows)
        wb_line = hlines!(ax, [wb.tHalf], color = :black, linestyle = (:dash, :dense), linewidth = 3.5)
        push!(handles, wb_line)
        push!(labels, world_best_label_tex(wb.project, signal))
    end

    hidedecorations!(ax, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)
    hidedecorations!(ax_top, grid = true, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

    ax.yticks = (yt_vals, yt_labels)
    ax_top.yticks = (yt_vals, yt_labels)

    Legend(
        fig[1, 2],
        handles,
        labels,
        legend_title_tex(signal),
        patchsize = (55, 30, 40, 80),
        patchlabelgap = 12,
        tellwidth = true,
        width = 450,
        rowgap = 10
    )

    out_dir = joinpath(base_dir, "figs")
    mkpath(out_dir)
    out_file = joinpath(out_dir, "sensitivity_vs_time_cases_$(signal_filename(signal)).png")
    save(out_file, fig, px_per_unit = 4)
    return fig, out_file
end

for signal in unique(sens_df.signal)
    signal_rows = filter(:signal => ==(signal), sens_df)
    world_best_rows = filter(:signal => ==(signal), world_bests_df)

    if nrow(signal_rows) < 3
        @warn "Skipping signal with insufficient radon points" signal
        continue
    end

    _, out_path = make_signal_plot(signal_rows, world_best_rows, signal)
    @info "Saved plot" signal out_path
end

