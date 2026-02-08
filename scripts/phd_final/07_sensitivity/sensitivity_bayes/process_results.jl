using DrWatson
@quickactivate "SNSensitivityEstimate"

using DataFramesMeta, CSV, CairoMakie, PrettyTables, LaTeXStrings, StatsBase

resultsdir = scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results/nu0bb")

d_nu0bb_tag1 = DataFrame()
for f in readdir(resultsdir)
    if occursin("tag=1_signal=bb0nu_foil_bulk", f)
        @show joinpath(resultsdir, f)
        d_nu0bb_tag1 = vcat(d_nu0bb_tag1, CSV.File(joinpath(resultsdir, f)) |> DataFrame)
    end
end
d_nu0bb_tag1 = d_nu0bb_tag1[1:10_000, :]


d_nu0bb_tag2 = DataFrame()
for f in readdir(resultsdir)
    if occursin("tag=2_signal=bb0nu_foil_bulk", f)
        @show joinpath(resultsdir, f)
        d_nu0bb_tag2 = vcat(d_nu0bb_tag2, CSV.File(joinpath(resultsdir, f)) |> DataFrame)
    end
end
d_nu0bb_tag2 = d_nu0bb_tag2[1:10_000, :]
    
d_nu0bbM1_tag1 = DataFrame()
for f in readdir(resultsdir)
    if occursin("tag=1_signal=bb0nuM1_foil_bulk", f)
        @show joinpath(resultsdir, f)
        d_nu0bbM1_tag1 = vcat(d_nu0bbM1_tag1, CSV.File(joinpath(resultsdir, f)) |> DataFrame)
    end
end
d_nu0bbM1_tag1 = d_nu0bbM1_tag1[1:10_000, :]

d_nu0bbM2_tag1 = DataFrame()
for f in readdir(resultsdir)
    if occursin("tag=1_signal=bb0nuM2_foil_bulk", f)
        @show joinpath(resultsdir, f)
        d_nu0bbM2_tag1 = vcat(d_nu0bbM2_tag1, CSV.File(joinpath(resultsdir, f)) |> DataFrame)
    end
end
d_nu0bbM2_tag1 = d_nu0bbM2_tag1[1:10_000, :]

d_RH050_tag1 = DataFrame()
for f in readdir(resultsdir)
    if occursin("tag=1_signal=RH050_foil_bulk", f)
        @show joinpath(resultsdir, f)
        d_RH050_tag1 = vcat(d_RH050_tag1, CSV.File(joinpath(resultsdir, f)) |> DataFrame)
    end
end



d_total = DataFrame(
    process = [L"0\nu\beta\beta", L"0\nu\beta\beta", L"0\nu\beta\beta M1", L"0\nu\beta\beta M2", L"$\nu_R\nu_L\beta\beta \, (K=0.5)$"],
    radon_level = [L"150\,$\mu$Bq/m$^3$", L"2\,$m$Bq/m$^3$", L"150\,$\mu$Bq/m$^3$", L"150\,$\mu$Bq/m$^3$", L"150\,$\mu$Bq/m$^3$"],
    sensitivity = [median(d_nu0bb_tag1.thalf), median(d_nu0bb_tag2.thalf), median(d_nu0bbM1_tag1.thalf), median(d_nu0bbM2_tag1.thalf), median(d_RH050_tag1.thalf)],
)

colors = Makie.wong_colors()
with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"$0\nu\beta\beta$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, nothing, 0, nothing))
    stephist!(a, d_nu0bb_tag1.thalf, bins = 100, color = colors[1], linewidth = 4)
    vlines!(a, [median(d_nu0bb_tag1.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median sensitivity \n %$(round(median(d_nu0bb_tag1.thalf), sigdigits = 3))yr $$")
    axislegend(a, position = :rt, patchsize = (45, 35))
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/bb0nu_tag_1_median_sensitivity.png"), f, px_per_unit = 2)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"$0\nu\beta\beta$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, nothing, 0, nothing))
    stephist!(a, d_nu0bb_tag2.thalf, bins = 100, color = colors[2], linewidth = 4)
    vlines!(a, [median(d_nu0bb_tag2.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median sensitivity \n %$(round(median(d_nu0bb_tag2.thalf), sigdigits = 3))yr $$")
    axislegend(a, position = :rt, patchsize = (45, 35))
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/bb0nu_tag_2_median_sensitivity.png"), f, px_per_unit = 2)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"Bayesian sensitivity to $0\nu\beta\beta$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, maximum(vcat(d_nu0bb_tag1.thalf, d_nu0bb_tag2.thalf))*1.3, 0, nothing))
    stephist!(a, d_nu0bb_tag1.thalf, bins = 100, color = colors[1], linewidth = 4)
    vlines!(a, [median(d_nu0bb_tag1.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median $T^{0\nu} \geq %$(round(median(d_nu0bb_tag1.thalf), sigdigits = 3)) \, \mathrm{yr}$ \n $(a_{Rn}=150 \, \mathrm{\mu Bq/m^3} )$")
    stephist!(a, d_nu0bb_tag2.thalf, bins = 100, color = colors[2], linewidth = 4)
    vlines!(a, [median(d_nu0bb_tag2.thalf)], color = :black, linewidth = 4, linestyle = :dash, label = L"median $T^{0\nu} \geq %$(round(median(d_nu0bb_tag2.thalf), sigdigits = 3)) \, \mathrm{yr}$ \n $(a_{Rn}=2 \, \mathrm{mBq/m^3} )$")
    axislegend(a, position = :rt, patchsize = (45, 35))

    
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/nu0bb_sensitivities_together.png"), f, px_per_unit = 2)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"Bayesian sensitivity to $0\nu\beta\beta\chi^0$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, nothing, 0, nothing))
    stephist!(a, d_nu0bbM1_tag1.thalf, bins = 100, color = colors[3], linewidth = 4)
    vlines!(a, [median(d_nu0bbM1_tag1.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median sensitivity \n %$(round(median(d_nu0bbM1_tag1.thalf), sigdigits = 3))yr $$")
    axislegend(a, position = :rt, patchsize = (45, 35))
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/bb0nuM1_tag_1_median_sensitivity.png"), f, px_per_unit = 2)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"Bayesian sensitivity to $0\nu\beta\beta\chi^0\chi^0$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, nothing, 0, nothing))
    stephist!(a, d_nu0bbM2_tag1.thalf, bins = 100, color = colors[4], linewidth = 4)
    vlines!(a, [median(d_nu0bbM2_tag1.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median sensitivity \n %$(round(median(d_nu0bbM2_tag1.thalf), sigdigits = 3))yr $$")
    axislegend(a, position = :rt, patchsize = (45, 35))
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/bb0nuM2_tag_1_median_sensitivity.png"), f, px_per_unit = 2)
    f
end

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 800), fontsize = 34,)
    a = Axis(f[1, 1], title = L"Bayesian sensitivity to $\nu_R\nu_L\beta\beta \, (K=0.5)$", xlabel = L"90% C.I. sensitivity (yr)$$", ylabel = L"number of pseudo-experiments $$", limits=(nothing, nothing, 0, nothing))
    stephist!(a, d_RH050_tag1.thalf, bins = 100, color = colors[5], linewidth = 4)
    vlines!(a, [median(d_RH050_tag1.thalf)], color = :red, linewidth = 4, linestyle = :dash, label = L"median sensitivity \n %$(round(median(d_RH050_tag1.thalf), sigdigits = 3))yr $$")
    axislegend(a, position = :rt, patchsize = (45, 35))
    save(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/figs/RH050_tag_1_median_sensitivity.png"), f, px_per_unit = 2)
    f
end

begin
    open(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results_summary.tex"), "w") do io
        pretty_table(
            io, 
            d_total, 
            backend = Val(:latex),
            header = ["Process", "Radon level", "Median sensitivity (yr)"],
        )
    end

    open(scriptsdir("phd_final/07_sensitivity/sensitivity_bayes/results_summary.md"), "w") do io
        pretty_table(
            io, 
            d_total, 
            backend = Val(:markdown),
            header = ["Process", "Radon level", "Median sensitivity (yr)"],
        )
    end
end
