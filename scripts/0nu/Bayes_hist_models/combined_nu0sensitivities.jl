using CSV, DataFramesMeta, CairoMakie, StatsBase

t = Float64[]
signal = "bb0nu_foil_bulk"

for f in readdir("data_perBkg_10keV_binning")
    if occursin(signal, f)
        d = CSV.File(joinpath("data_perBkg_10keV_binning", f)) |> DataFrame
        append!(t, d.thalf)
    end
end

let
    f = Figure()
    a = Axis(f[1,1], xlabel = "sensitivity (yr)", ylabel = ("n pseudo-experiments"), limits = (nothing, nothing, 0, nothing))
    p1 = stephist!(a,t, bins = 40, label = "data")
    p2 = vlines!(
        a, 
        [mean(t)], 
        color = Makie.wong_colors()[2],
        linewidth = 3
    )
    p3 = vlines!(
        a, 
        [median(t)], 
        color = Makie.wong_colors()[3],
        linewidth = 3
    )

    Legend(f[1,2], [p1, p2, p3], ["data", "mean = $(round(mean(t), sigdigits = 3)) yr", "median = $(round(median(t), sigdigits = 3)) yr"])
    save("fig.png", f)
    f
end