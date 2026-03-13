using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using UnROOT

f = ROOTFile("data/data/final_phd/comparison_E_corrections/combined_p3.root")

df = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame

@chain df begin
    @subset! :E_f1 .> 0
    @subset! :E_f2 .> 0
    @subset! :E_f_optical1 .> 0
    @subset! :E_f_optical2 .> 0
    @subset! :E_f_optical_loss1 .> 0
    @subset! :E_f_optical_loss2 .> 0
end

df_e = DataFrame(
    E_f = vcat(df.E_f1, df.E_f2),
    E_f_optical = vcat(df.E_f_optical1, df.E_f_optical2),
    E_f_optical_loss = vcat(df.E_f_optical_loss1, df.E_f_optical_loss2)
)

binning = 0:50:5000

h_E_f = Hist1D(df_e.E_f; binedges=binning)
h_E_f_optical = Hist1D(df_e.E_f_optical; binedges=binning)
h_E_f_optical_loss = Hist1D(df_e.E_f_optical_loss; binedges=binning)

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 600), fontsize = 32)
    a = Axis(f[1, 1], title = "Comparison of various energy corrections; phase 3", xlabel = "Single-electron energy (keV)", ylabel = "Counts / 50 keV")
    stephist!(a, h_E_f, label = L"$E_f$ final energy", linewidth = 3)
    stephist!(a, h_E_f_optical, label = L"$E_D$ optical corrections", linewidth = 3)
    stephist!(a, h_E_f_optical_loss, label = L"$E_i$ initial energy", linewidth = 3)
    axislegend(a, position = :rt)
    save("scripts/phd_final/08_data/figs/compare_corrections/comparison_single_E_corrections.png", f, px_per_unit = 4)
    f
end

df_e_sum = DataFrame(
    E_f_sum = df.E_f1 + df.E_f2,
    E_f_optical_sum = df.E_f_optical1 + df.E_f_optical2,
    E_f_optical_loss_sum = df.E_f_optical_loss1 + df.E_f_optical_loss2
)

h_E_f_sum = Hist1D(df_e_sum.E_f_sum; binedges=binning)
h_E_f_optical_sum = Hist1D(df_e_sum.E_f_optical_sum; binedges=binning)
h_E_f_optical_loss_sum = Hist1D(df_e_sum.E_f_optical_loss_sum; binedges=binning)

with_theme(theme_latexfonts()) do
    f = Figure(size = (1200, 600), fontsize = 32)
    a = Axis(f[1, 1], title = "Comparison of various energy corrections; phase 3", xlabel = "Summed-electron energy (keV)", ylabel = "Counts / 50 keV")
    stephist!(a, h_E_f_sum, label = L"$E_f$ final energy", linewidth = 3)
    stephist!(a, h_E_f_optical_sum, label = L"$E_D$ optical corrections", linewidth = 3)
    stephist!(a, h_E_f_optical_loss_sum, label = L"$E_i$ initial energy", linewidth = 3)
    a.xticks = 0:500:5000
    axislegend(a, position = :rt)
    save("scripts/phd_final/08_data/figs/compare_corrections/comparison_two_E_corrections.png", f, px_per_unit = 4)
    f
end

function get_fwhm(h::Hist1D)
    bin_centers = bincenters(h)
    max_bin_id = argmax(bincounts(h))
    half_max = bincounts(h)[max_bin_id] / 2
    left_idx = findlast(x -> x < half_max, bincounts(h)[1:max_bin_id])
    right_idx = findfirst(x -> x < half_max, bincounts(h)[max_bin_id:end]) + max_bin_id - 1
    fwhm = bin_centers[right_idx] - bin_centers[left_idx]
    return fwhm
end


function fwhm(h::Hist1D)
    counts = bincounts(h)
    edges = binedges(h)[1]

    # bin centers
    centers = bincenters(h)

    # maximum
    imax = argmax(counts)
    halfmax = counts[imax] / 2

    # ---- find left crossing ----
    left = nothing
    for i in reverse(1:imax)
        if counts[i] < halfmax
            # linear interpolation
            x1, x2 = centers[i], centers[i+1]
            y1, y2 = counts[i], counts[i+1]
            left = x1 + (halfmax - y1) * (x2 - x1) / (y2 - y1)
            break
        end
    end

    # ---- find right crossing ----
    right = nothing
    for i in imax:length(counts)-1
        if counts[i+1] < halfmax
            x1, x2 = centers[i], centers[i+1]
            y1, y2 = counts[i], counts[i+1]
            right = x1 + (halfmax - y1) * (x2 - x1) / (y2 - y1)
            break
        end
    end

    if left === nothing || right === nothing
        error("Could not determine FWHM (histogram may not cross half max).")
    end

    return right - left
end

fwhm_E_f_sum = fwhm(h_E_f_sum)
fwhm_E_f_optical_sum = fwhm(h_E_f_optical_sum)
fwhm_E_f_optical_loss_sum = fwhm(h_E_f_optical_loss_sum)

