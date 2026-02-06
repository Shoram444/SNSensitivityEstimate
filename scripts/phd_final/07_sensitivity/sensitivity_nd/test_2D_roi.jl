using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist

println("loaded pkgs")

α= 1.64485362695147


d_s = CSV.File("scripts/phd_final/07_sensitivity/sensitivity_nd/figs/for_miro_ND_signal_2Dhist_sumE_vs_phi.csv") |> DataFrame
d_b = CSV.File("scripts/phd_final/07_sensitivity/sensitivity_nd/figs/for_miro_ND_background_2Dhist_sumE_vs_phi.csv") |> DataFrame

d_joined = innerjoin(d_s, d_b, on=[:Esum, :phi])

@transform! d_joined :s_b = :eff ./ get_FC.(:bkg_counts, α; approximate = "table")
@transform! d_joined :b_excl = get_FC.(:s_b, α; approximate = "table")

sort!(d_joined, :s_b, rev=true)




function accumulate_sens(df::DataFrame)

    e1 = df.eff[1]
    b1 = df.b_excl[1]
    s1 = get_tHalf(SNparams, e1, b1, α; approximate="table")

    selected = [1]

    for i in 2:nrow(df)
        etot = e1 + df.eff[i]
        btot = b1 + df.b_excl[i]

        sens = get_tHalf(SNparams, etot, btot, α; approximate="table")

        if sens > s1
            
            e1 = etot
            b1 = btot
            s1 = sens
            @show e1, b1, s1
            push!(selected, i)
        end
    end

    return (e1=e1, b1=b1, sens=s1, selected_bins=selected)
end

res = accumulate_sens(d_joined)


function opt1(df::DataFrame)

    e1 = df.eff[1]
    b1 = df.b_excl[1]

    s1 = get_tHalf(SNparams, e1, b1, α; approximate="table")

    selected = [1]

    for i in 2:nrow(df)
        etot = e1 + df.eff[i]
        btot = b1 + df.b_excl[i]

        snew = get_tHalf(SNparams, etot, btot, α; approximate="table")

        Δs = snew - s1

        if Δs > 0
            e1 = etot
            b1 = btot
            s1 = snew
            push!(selected, i)

            # @show i, e1, b1, s1
        end
    end

    return (
        e1 = e1,
        b1 = b1,
        sens = s1,
        selected_bins = selected
    )
end

res_opt1 = opt1(d_joined)

function opt2(df::DataFrame)

    e1 = df.eff[1]
    b1 = df.b_excl[1]

    s1 = get_tHalf(SNparams, e1, b1, α; approximate="table")
    r1 = s1 / e1

    selected = [1]

    for i in 2:nrow(df)
        e2 = df.eff[i]
        b2 = df.b_excl[i]

        etot = e1 + e2
        btot = b1 + b2

        snew = get_tHalf(SNparams, etot, btot, α; approximate="table")
        rnew = snew / etot

        if rnew > r1
            e1 = etot
            b1 = btot
            s1 = snew
            r1 = rnew
            push!(selected, i)

            # @show i, e1, b1, s1, r1
        else
            break
        end
    end

    return (
        e1 = e1,
        b1 = b1,
        sens = s1,
        selected_bins = selected
    )
end

res_opt2 = opt2(d_joined)



function pareto_prune(df::DataFrame)

    # Sort by increasing background, decreasing efficiency
    sdf = sort(df, [:b_excl], rev = false)

    keep = BitVector(undef, nrow(sdf))
    best_eff = -Inf

    for i in 1:nrow(sdf)
        if sdf.eff[i] > best_eff
            keep[i] = true
            best_eff = sdf.eff[i]
        else
            keep[i] = false
        end
    end

    return sdf[keep, :]
end

res_pareto = pareto_prune(d_joined)


function best_subset_beam(
    df::DataFrame,
    SNparams,
    α;
    K = 100,
    approximate = "table",
    max_bins = nrow(df)
)

    # State: (indices, eff, bkg, sens)
    states = [ (Int[], 0.0, 0.0, 0.0) ]

    for i in 1:min(max_bins, nrow(df))
        new_states = Vector{Tuple{Vector{Int},Float64,Float64,Float64}}()

        for (idx, eff, bkg, _) in states
            eff2 = eff + df.eff[i]
            bkg2 = bkg + df.b_excl[i]

            sens2 = get_tHalf(SNparams, eff2, bkg2, α; approximate)

            push!(new_states, (vcat(idx, i), eff2, bkg2, sens2))
        end

        # Also allow "skip bin"
        append!(new_states, states)

        # Keep top K
        sort!(new_states, by = x -> x[4], rev = true)
        states = new_states[1:min(K, length(new_states))]
    end

    best = argmax(x -> x[4], states)
    idx, eff, bkg, sens = states[best]

    return (
        sens = sens,
        eff  = eff,
        bkg  = bkg,
        bins = df[idx, :]
    )
end
res_beam = best_subset_beam(d_joined, SNparams, α; K=200, max_bins=50)