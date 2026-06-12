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

colors = ["#041E42", "#BE4D00", "#951272", "#006630", "#005C8A", "#FFB948", "#605643", "#302D23"]
line_styles = [:solid, :dash, :dashdot, :dot, :dashdotdot]

signal_names = Dict(
    "bb0nu_foil_bulk" => "0nu beta beta",
    "bb0nuM1_foil_bulk" => "0nu beta beta chi0",
    "bb0nuM2_foil_bulk" => "0nu beta beta chi0 chi0",
    "RH050_foil_bulk" => "nuL nuR beta beta",
    "RH040_foil_bulk" => "nuL nuR beta beta",
)

function signal_title(signal::AbstractString)
    return get(signal_names, signal, replace(signal, "_" => " "))
end

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

function radon_label(radon_level_bq_m3)
    radon_mBq_m3 = radon_level_bq_m3 * 1_000.0
    return @sprintf("Rn = %.2f mBq/m^3", radon_mBq_m3)
end

function radon_label_tex(radon_level_bq_m3)
    if isapprox(radon_level_bq_m3, 0.03821; atol = 1e-6)
        return L"$$\ Rn @ 38.21\ \mathrm{mBq/m^3}"
    elseif isapprox(radon_level_bq_m3, 0.00108; atol = 1e-6)
        return L"$$\ Rn @ 1.08\ \mathrm{mBq/m^3}"
    elseif isapprox(radon_level_bq_m3, 0.00015; atol = 1e-6)
        return L"$$\ Rn @ 0.15\ \mathrm{mBq/m^3}"
    end
    return L"$$\ Rn\ level"
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

function compute_sensitivity_curve(eff, bkg, times, t_base; efffactor = 1.0, extra_bkg = 0.0)
    α = 1.64485362695147
    t_halfs = Float64[]
    sizehint!(t_halfs, length(times))
    for t in times
        scaled_bkg = bkg / t_base * t
        t_half = get_tHalf(
            SNparams["W"],
            SNparams["foilMass"],
            SNparams["Nₐ"],
            t,
            SNparams["a"],
            eff * efffactor,
            scaled_bkg + extra_bkg,
            α;
            approximate = "table",
        )
        push!(t_halfs, t_half)
    end
    return t_halfs
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

function decade_ticks(y_low, y_high)
    p_low = Int(floor(log10(y_low)))
    p_high = Int(ceil(log10(y_high)))
    vals = [10.0^p for p in p_low:p_high]
    labels = [latexstring("10^{", string(p), "}") for p in p_low:p_high]
    return vals, labels
end

function signal_filename(signal::AbstractString)
    return replace(signal, r"[^A-Za-z0-9]+" => "_")
end

function make_signal_plot(signal_rows::DataFrame, world_best_rows::DataFrame, signal::AbstractString)
    sort!(signal_rows, :radon_level, rev = true)

    curves = [compute_sensitivity_curve(r.eff, r.bkgCounts, times, t_base) for r in eachrow(signal_rows)]
    y_low, y_high = y_limits_log(curves)
    yt_vals, yt_labels = decade_ticks(y_low, y_high)

    fig = Figure(
        size = (950*1.3, 600),
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

    handles = Any[]
    labels = Any[]

    for (idx, row) in enumerate(eachrow(signal_rows))
        color = colors[mod1(idx, length(colors))]
        style = line_styles[mod1(idx, length(line_styles))]
        p = lines!(ax, exposures, curves[idx], color = color, linewidth = 5.5, linestyle = style)
        lines!(ax_top, times, curves[idx], color = color, linewidth = 5.5, linestyle = style)
        push!(handles, p)
        push!(labels, radon_label_tex(row.radon_level))
    end

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
        patchlabelgap  = 20,
        tellwidth = true,
        width = 450
    )

    out_dir = joinpath(base_dir, "figs")
    mkpath(out_dir)
    out_file = joinpath(out_dir, "sensitivity_vs_time_$(signal_filename(signal)).png")
    save(out_file, fig, px_per_unit = 4)
    return fig, out_file
end

for signal in unique(sens_df.signal)
    signal_rows = filter(:signal => ==(signal), sens_df)
    world_best_rows = filter(:signal => ==(signal), world_bests_df)
    _, out_path = make_signal_plot(signal_rows, world_best_rows, signal)
    @info "Saved plot" signal out_path
end

