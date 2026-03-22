using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using StatsBase
using UnROOT, ColorSchemes, StatsPlots

in1 = CSV.File(datadir("Ovidiu_spectra_bsm/3DDistributions_NoSterile_82Se.csv")) |> DataFrame
in2 = CSV.File(datadir("Ovidiu_spectra_bsm/3DDistributions_Sterile_All_Masses_82Se.csv")) |> DataFrame


gdf_in2 = @chain in2 begin
    @subset! :e1_MeV .>= 0.01
    @subset! :e2_MeV .>= 0.01
    @rtransform :Gamma_sterile = ifelse(:Gamma_sterile < 1e-9, 0.0, :Gamma_sterile)
    @groupby :mN_MeV 
end

gdf_in1 = @chain in1 begin
    @transform :e1_MeV = round.(:e1_MeV, digits = 3)
    @transform :e2_MeV = round.(:e2_MeV, digits = 3)
    @transform :Gamma_SM_N0 = round.(:Gamma_SM_N0, digits = 3)
    @transform :Gamma_Majoron_m1n1 = round.(:Gamma_Majoron_m1n1, digits = 3)
    @transform :Gamma_Majoron_m2n3 = round.(:Gamma_Majoron_m2n3, digits = 3)
    @subset :e1_MeV .>= 0.01
    @subset :e2_MeV .>= 0.01
    @transform :esum = :e1_MeV + :e2_MeV
    @groupby :e1_MeV :e2_MeV
    @combine begin
        :Gamma_SM_N0 = sum(:Gamma_SM_N0)
        :Gamma_Majoron_m1n1 = sum(:Gamma_Majoron_m1n1)
        :Gamma_Majoron_m2n3 = sum(:Gamma_Majoron_m2n3)
    end
    # @orderby :esum
end

StatsPlots.plot(gdf_in1.e1_MeV, gdf_in1.e2_MeV, gdf_in1.Gamma_SM_N0)

nbins = length(unique(gdf_in1.e1_MeV))
h2d = Hist2D(;binedges = (range(0.01, 2.98, length = nbins + 1), range(0.01, 2.98, length = nbins + 1)))
for row in eachrow(gdf_in1)
    push!(h2d, row.e1_MeV, row.e2_MeV, row.Gamma_SM_N0)
end




"""
    project_2d_probability_to_esum(df::DataFrame, e1_col::Symbol, e2_col::Symbol, rate_col::Symbol)

Project a 2D probability distribution P(e1, e2) to P(esum) analytically.
For each esum, sums rates from all (e1, e2) pairs where e1 + e2 = esum.

Returns a tuple of (esum_values, probabilities).
"""
function project_2d_probability_to_esum(
    df::DataFrame,
    e1_col::Symbol,
    e2_col::Symbol,
    rate_col::Symbol;
    esum_digits::Union{Nothing,Int}=3,
    esum_binwidth::Union{Nothing,Float64}=nothing,
)
    esum = Float64.(df[!, e1_col] .+ df[!, e2_col])

    # Stabilize floating-point artifacts: either explicit binning or decimal quantization.
    if !isnothing(esum_binwidth)
        bw = esum_binwidth
        emin = minimum(esum)
        esum = emin .+ bw .* round.((esum .- emin) ./ bw)
    elseif !isnothing(esum_digits)
        esum = round.(esum, digits=esum_digits)
    end

    tmp = DataFrame(esum=esum, rate=Float64.(df[!, rate_col]))

    proj = @chain tmp begin
        @groupby :esum
        @combine :rate = sum(:rate)
        @orderby :esum
    end

    return proj.esum, proj.rate
end


"""
    sample_esum_from_2d_probability(df::DataFrame, e1_col::Symbol, e2_col::Symbol, rate_col::Symbol;
                                    nsamples::Int=200_000, normalize::Bool=true)

Draw samples from discrete weighted points (e1, e2, rate), then return sampled esum = e1 + e2.
If `normalize=true`, rates are normalized internally to probabilities.
"""
function sample_esum_from_2d_probability(
    df::DataFrame,
    e1_col::Symbol,
    e2_col::Symbol,
    rate_col::Symbol;
    nsamples::Int=200_000,
    normalize::Bool=true,
    jitter_within_bin::Bool=true,
    e1_binwidth::Union{Nothing,Float64}=nothing,
    e2_binwidth::Union{Nothing,Float64}=nothing,
)
    w = Float64.(df[!, rate_col])
    if any(<(0.0), w)
        error("rate_col contains negative values, cannot sample probabilities")
    end

    if normalize
        s = sum(w)
        s > 0 || error("sum of rates is zero, cannot sample")
        w ./= s
    end

    idx = StatsBase.sample(1:nrow(df), Weights(w), nsamples)

    e1 = Float64.(df[idx, e1_col])
    e2 = Float64.(df[idx, e2_col])

    if jitter_within_bin
        if isnothing(e1_binwidth)
            u1 = sort(unique(Float64.(df[!, e1_col])))
            if length(u1) > 1
                e1_binwidth = median(diff(u1))
            end
        end
        if isnothing(e2_binwidth)
            u2 = sort(unique(Float64.(df[!, e2_col])))
            if length(u2) > 1
                e2_binwidth = median(diff(u2))
            end
        end

        if !isnothing(e1_binwidth)
            e1 .+= (rand(nsamples) .- 0.5) .* e1_binwidth
        end
        if !isnothing(e2_binwidth)
            e2 .+= (rand(nsamples) .- 0.5) .* e2_binwidth
        end
    end

    return e1 .+ e2
end

esum_vals, rate_SM = project_2d_probability_to_esum(
    in1,
    :e1_MeV,
    :e2_MeV,
    :Gamma_SM_N0;
    esum_binwidth=0.01,
)

_, rate_Majoron_m1n1 = project_2d_probability_to_esum(
    in1,
    :e1_MeV,
    :e2_MeV,
    :Gamma_Majoron_m1n1;
    esum_binwidth=0.01,
)

_, rate_Majoron_m2n3 = project_2d_probability_to_esum(
    in1,
    :e1_MeV,
    :e2_MeV,
    :Gamma_Majoron_m2n3;
    esum_binwidth=0.01,
)


colors = ["#003865", "#FFB948", "#52473B", "#9A3A06", "#951272", "#006630", "#005C8A", "#FFB948"]
with_theme(theme_latexfonts()) do  
    rate_SM ./= maximum(rate_SM)
    rate_Majoron_m1n1 ./= maximum(rate_Majoron_m1n1)
    rate_Majoron_m2n3 ./= maximum(rate_Majoron_m2n3)
    f = Figure(size = (1200, 800), fontsize = 28, figure_padding = 35)
    a = Axis(f[1,1], title = "Summed-electron energy distributions of DBDs with Majoron emission", xlabel = L"$E_{sum}$ (MeV)", ylabel = L"rate (a.u.) $$", limits = (0, 4, 0, 1.05))

    CairoMakie.lines!(a, esum_vals, rate_SM, label = L"\text{SM:}~ 2\nu\beta\beta", color = colors[1], linewidth = 4)
    CairoMakie.band!(a, esum_vals, zero.(rate_SM), rate_SM, color = (colors[1], 0.15))
    CairoMakie.lines!(a, esum_vals, rate_Majoron_m1n1, label = L"0\nu\beta\beta\chi^0 ~(n=1)", color = colors[2], linewidth = 4, linestyle = :dash)
    CairoMakie.lines!(a, esum_vals, rate_Majoron_m2n3, label = L"0\nu\beta\beta\chi^0\chi^0 ~(n=3)", color = colors[3], linewidth = 4, linestyle = :dashdot)
    
    axislegend(a, position = :rt, patchsize = (55, 30, 40, 80), patchlabelgap  = 20)
    save(scriptsdir("phd_final/02_BSM_theo/figs/esum_majoron.png"), f, px_per_unit = 3)
    f
end

## 

esum_05, sterile_05 = project_2d_probability_to_esum(
    DataFrame(gdf_in2[1]),
    :e1_MeV,
    :e2_MeV,
    :Gamma_sterile;
    esum_binwidth=0.01,
)

esum_10, sterile_10 = project_2d_probability_to_esum(
    DataFrame(gdf_in2[6]),
    :e1_MeV,
    :e2_MeV,
    :Gamma_sterile;
    esum_binwidth=0.01,
)

esum_15, sterile_15 = project_2d_probability_to_esum(
    DataFrame(gdf_in2[11]),
    :e1_MeV,
    :e2_MeV,
    :Gamma_sterile;
    esum_binwidth=0.01,
)

with_theme(theme_latexfonts()) do  
    sterile_05 ./= maximum(sterile_05)
    sterile_10 ./= maximum(sterile_10)
    sterile_15 ./= maximum(sterile_15)
    f = Figure(size = (1200, 800), fontsize = 28, figure_padding = 35)
    a = Axis(f[1,1], title = "Summed-electron energy distributions of DBDs with sterile neutrino", xlabel = L"$E_{sum}$ (MeV)", ylabel = L"rate (a.u.) $$", limits = (0, 4, 0, 1.05))

    CairoMakie.lines!(a, esum_vals, rate_SM, label = L"\text{SM:}~ 2\nu\beta\beta", color = colors[1], linewidth = 4)
    CairoMakie.band!(a, esum_vals, zero.(rate_SM), rate_SM, color = (colors[1], 0.15))
    CairoMakie.lines!(a, esum_05, sterile_05, label = L"N\nu\beta\beta ~(m_N=0.5~\text{MeV})", color = colors[4], linewidth = 4, linestyle = :dash)
    CairoMakie.lines!(a, esum_10, sterile_10, label = L"N\nu\beta\beta ~(m_N=1.0~\text{MeV})", color = colors[5], linewidth = 4, linestyle = :dashdot)
    CairoMakie.lines!(a, esum_15, sterile_15, label = L"N\nu\beta\beta ~(m_N=1.5~\text{MeV})", color = colors[6], linewidth = 4, linestyle = :dot)

    axislegend(a, position = :rt, patchsize = (55, 30, 40, 80), patchlabelgap  = 20)
    save(scriptsdir("phd_final/02_BSM_theo/figs/esum_sterile.png"), f, px_per_unit = 3)
    f
end


for g in gdf_in2
    m_n = unique(g.mN_MeV)[1]
    CSV.write(datadir("Ovidiu_spectra_bsm/sterile_$m_n.csv"), DataFrame(e1_MeV = g.e1_MeV, e2_MeV = g.e2_MeV, Gamma_sterile = g.Gamma_sterile), header =false)
end

