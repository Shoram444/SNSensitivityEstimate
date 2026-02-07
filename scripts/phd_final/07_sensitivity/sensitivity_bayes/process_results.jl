using DrWatson
@quickactivate "SNSensitivityEstimate"

using DataFramesMeta, CSV, CairoMakie, PrettyTables, LaTeXStrings, StatsBase

d_nu0bb_tag1 = CSV.File(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/nu0bb/nu0bb_tag1_small_run.csv")) |> DataFrame

d_total = DataFrame(
    process = [L"0\nu\beta\beta"],
    radon_level = [L"150\,$\mu$Bq/m$^3$"],
    sensitivity = [median(d_nu0bb_tag1.thalf)],
)

colors = Makie.wong_colors()
with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"$0\nu\beta\beta$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, nothing, 0, nothing))
    stephist!(a, d_nu0bb_tag1.thalf, bins = 20, color = colors[1], linewidth = 4)
    vlines!(a, [median(d_nu0bb_tag1.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median sensitivity $$")
    axislegend(a, position = :rt, patchsize = (45, 35))
    f
end