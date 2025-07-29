using CairoMakie, CSV, DataFramesMeta, UnROOT, LaTeXStrings, FHist, StatsBase, KernelDensity
colors          = ["#003865", "#FFB948", "#52473B", "#9A3A06", "#5AB4C0"]

f_m1 = ROOTFile(joinpath("scripts/Medex25/theo_data/m1.root"))
df_m1 = LazyTree(f_m1, "tree", keys(f_m1["tree"])) |> DataFrame
# @transform! df_m1 :esum = :simuEnergy1 .+ :simuEnergy2
# @subset! df_m1 :esum .>= 10.0

f_m2 = ROOTFile("scripts/Medex25/theo_data/m2.root")
df_m2 = LazyTree(f_m2, "tree", keys(f_m2["tree"])) |> DataFrame


f_rh = ROOTFile("scripts/Medex25/theo_data/rh.root")
df_rh = LazyTree(f_rh, "tree", keys(f_rh["tree"])) |> DataFrame

f_nu0 = ROOTFile("scripts/Medex25/theo_data/nu0.root")
df_nu0 = LazyTree(f_nu0, "tree", keys(f_nu0["tree"])) |> DataFrame


f_bb = ROOTFile("scripts/Medex25/theo_data/bb.root")
df_bb = LazyTree(f_bb, "tree", keys(f_bb["tree"])) |> DataFrame

bins = range(0, 1, 50)
bc = midpoints(bins) |> collect
bc[end] = 1
bc[1] = 0

hm1 = Hist1D(df_m1.esum ./ 3000 ; binedges = bins) |> normalize
hm2 = Hist1D(df_m2.esum ./ 3000 ; binedges = bins) |> normalize
hrh = Hist1D(df_rh.esum ./ 3000 ; binedges = bins) |> normalize
hbb = Hist1D(df_bb.esum ./ 3000 ; binedges = bins) |> normalize
hnu0 = Hist1D(df_nu0.esum ./ 3000 ; binedges = bins) |> normalize

bincountsm1 = bincounts(hm1) ./ maximum(bincounts(hm1))
bincountsm2 = bincounts(hm2) ./ maximum(bincounts(hm2))
bincountsrh = bincounts(hrh) ./ maximum(bincounts(hrh))
bincountsbb = bincounts(hbb) ./ maximum(bincounts(hbb))
bincountsnu0 = bincounts(hnu0) ./ maximum(bincounts(hnu0))
bincountsm1[end] = 0
bincountsm2[end] = 0
bincountsrh[end] = 0
bincountsbb[end] = 0
bincountsnu0[end] = 0

## bb vs nu0
with_theme(theme_latexfonts()) do
    f = Figure(
        size = (1700, 600), 
        fontsize = 42,  
        figure_padding = 45,
    )
    a = Axis(
        f[1,1], 
        xlabel = L"summed 2-electron energy $E_1 + E_2$",
        ylabel = L"rate (a.u.) $$",
        title = L"theoretical spectra of double beta decay processes $$",
        limits = (0,1.08,0,1.05),
        xticks = (0:1, [L"0", L"Q"]),
        yticks = (0:1, [L"0", L"1"]),
    )

    p1 = lines!(a, [1, 1], [0, 0.3], label = L" $0\nu\beta\beta$", linewidth = 8, color = :black)
    p2 = lines!(a, bc, bincountsbb, label = L" $2\nu\beta\beta$", linewidth = 8, color = colors[2], linestyle = (:dot, :dense))
    # lines!(a, bc, bincountsm2, label = L" $0\nu\beta\beta\chi^2\chi^2$, n=3", linewidth = 8, color = colors[4])
    
    Legend(
        f[1,2], 
        [p1,p2], 
        [L" $0\nu\beta\beta$             ", L" $2\nu\beta\beta$ "], 
        patchsize = (80, 50),
        valign = :top,
        halign = :left,
        width = 400,
        # patchlabelgap = 80
    )
    # xticks!(a, (0,1), (0, "Q"))
    # axislegend(a, patchsize = (40, 50), position = :cb, nbanks = 2)
    hidedecorations!(a, label = false,  ticklabels = false, ticks = false,)
    save("scripts/Medex25/figs/theo/bb_vs_nu0.png",f, px_per_unit = 5)
    f
end

## all
with_theme(theme_latexfonts()) do
    f = Figure(
        size = (1700, 600), 
        fontsize = 42,  
        figure_padding = 45,
    )
    a = Axis(
        f[1,1], 
        xlabel = L"summed 2-electron energy $E_1 + E_2$ ",
        ylabel = L"rate (a.u.) $$",
        title = L"theoretical spectra of double beta decay processes $$",
        limits = (0,1.08,0,1.05),
        xticks = (0:1, [L"0", L"Q"]),
        yticks = (0:1, [L"0", L"1"]),
    )

    p1 = lines!(a, [1, 1], [0, 0.3], label = L" $0\nu\beta\beta$", linewidth = 8, color = :black)
    p2 = lines!(a, bc, bincountsbb, label = L" $2\nu\beta\beta$", linewidth = 8, color = colors[2], linestyle = (:dot, :dense))
    p3 = lines!(a, bc, bincountsm1, label = L" $0\nu\beta\beta\chi^0, n=1$", linewidth = 8, color = colors[4], linestyle = (:dash, :dense))
    p4 = lines!(a, bc, bincountsm2, label = L" $0\nu\beta\beta\chi^0\chi^0, n=3$", linewidth = 8, color = colors[1], linestyle = (:dashdot, :dense))
    # lines!(a, bc, bincountsm2, label = L" $0\nu\beta\beta\chi^2\chi^2$, n=3", linewidth = 8, color = colors[4])
    
    Legend(
        f[1,2], 
        [p1,p2, p3, p4], 
        [L" $0\nu\beta\beta$ ", L" $2\nu\beta\beta$ ", L" $0\nu\beta\beta\chi^0, n=1$ ", L" $0\nu\beta\beta\chi^0\chi^0, n=3$"], 
        patchsize = (80, 50),
        valign = :top,
        width = 400
    )
    # axislegend(a, patchsize = (40, 50), position = :cb, nbanks = 2)
    hidedecorations!(a, label = false,  ticklabels = false, ticks = false,)
    save("scripts/Medex25/figs/theo/esum.png",f, px_per_unit = 5)
    f
end


### single
bins = range(0, 1, 80)
bc = midpoints(bins) |> collect
bc[end] = 1
bc[1] = 0

bbnu0 = CSV.File("scripts/Medex25/theo_data/bbnu0.csv", delim = ",", ) |> DataFrame
bbnu0 = bbnu0[1:end-1, :]
@transform! bbnu0 :e = :e ./ maximum(:e)
@transform! bbnu0 :dg = :dg ./ maximum(:dg)

hm1Single = Hist1D(df_m1.simuEnergy1 ./ 3000 ; binedges = bins) |> normalize
hm2Single = Hist1D(df_m2.simuEnergy1 ./ 3000 ; binedges = bins) |> normalize
hrhSingle = Hist1D(df_rh.simuEnergy1 ./ 3000 ; binedges = bins) |> normalize
hnu0Single = Hist1D(df_nu0.simuEnergy1 ./ 3000 ; binedges = bins) |> normalize
bb2Single = Hist1D(df_bb.simuEnergy1 ./ 3000 ; binedges = bins) |> normalize

bincountsm1Single = bincounts(hm1Single) ./ maximum(bincounts(hm1Single))
bincountsm2Single = bincounts(hm2Single) ./ maximum(bincounts(hm2Single))
bincountsrhSingle = bincounts(hrhSingle) ./ maximum(bincounts(hrhSingle))
bincountsbbSingle = bincounts(bb2Single) ./ maximum(bincounts(bb2Single))

bincountsnu0Single = bincounts(hnu0Single) ./ maximum(bincounts(hnu0Single))
bincountsnu0Single = bincountsnu0Single .- minimum(bincountsnu0Single)
bincountsnu0Single = bincountsnu0Single ./ maximum(bincountsnu0Single)


with_theme(theme_latexfonts()) do
    f = Figure(
        size = (1700, 600), 
        fontsize = 42,  
        figure_padding = 45,
    )
    a = Axis(
        f[1,1], 
        xlabel = L"single electron energy $E_i$ ",
        ylabel = L"rate (a.u.) $$",
        title = L"theoretical spectra of double beta decay processes $$",
        limits = (0,1.08,0,1.05),
        xticks = (0:1, [L"0", L"Q"]),
        yticks = (0:1, [L"0", L"1"]),
    )

    p1 = lines!(a, bbnu0.e, bbnu0.dg, label = L" $0\nu\beta\beta$", linewidth = 8, color = :black)
    p2 = lines!(a, bc, bincountsbbSingle, label = L" $2\nu\beta\beta$", linewidth = 8, color = colors[2], linestyle = (:dot, :dense))
    p3 = lines!(a, bc, bincountsm1Single, label = L" $0\nu\beta\beta\chi^0, n=1$", linewidth = 8, color = colors[4], linestyle = (:dash, :dense))
    p4 = lines!(a, bc, bincountsm2Single, label = L" $0\nu\beta\beta\chi^0\chi^0, n=3$", linewidth = 8, color = colors[1], linestyle = (:dashdot, :dense))
    p5 = lines!(a, bc, bincountsrhSingle, label = L" $0\nu\beta\beta$ RH", linewidth = 8, color = colors[5], linestyle = (:dashdotdot, :dense))
    # lines!(a, bc, bincountsm2, label = L" $0\nu\beta\beta\chi^2\chi^2$, n=3", linewidth = 8, color = colors[4])
    
    Legend(
        f[1,2], 
        [p1,p2, p3, p4, p5], 
        [L" $0\nu\beta\beta$ ", L" $2\nu\beta\beta$ ", L" $0\nu\beta\beta\chi^0, n=1$ ", L" $0\nu\beta\beta\chi^0\chi^0, n=3$", L" $0\nu\beta\beta$ RH"], 
        patchsize = (80, 50),
        valign = :top,
        width = 400
    )
    # axislegend(a, patchsize = (40, 50), position = :cb, nbanks = 2)
    hidedecorations!(a, label = false,  ticklabels = false, ticks = false,)
    save("scripts/Medex25/figs/theo/esingle.png",f, px_per_unit = 5)
    f
end


### phi
using MPGenbb
thetas_RH = [rad2deg(sample_theta_dif(0.4)) for _=1:1e8]
h_RH = Hist1D(thetas_RH; binedges=0:5:180) 
bc_RH = midpoints(binedges(h_RH)) |> collect
ys_RH = bincounts(h_RH) ./ maximum(bincounts(h_RH))


binsPhi = range(0, 180, 50)
bcPhi = midpoints(binsPhi) |> collect
bcPhi[end] = 180
bcPhi[1] = 0

hm1Phi = Hist1D(df_m1.theta ; binedges = binsPhi) |> normalize
hm2Phi = Hist1D(df_m2.theta ; binedges = binsPhi) |> normalize
hrhPhi = Hist1D(df_rh.theta ; binedges = binsPhi) |> normalize
hbbPhi = Hist1D(df_bb.theta ; binedges = binsPhi) |> normalize
hnu0Phi = Hist1D(df_nu0.theta ; binedges = binsPhi) |> normalize

bincountsm1Phi = bincounts(hm1Phi) ./ maximum(bincounts(hm1Phi))
bincountsm2Phi = bincounts(hm2Phi) ./ maximum(bincounts(hm2Phi))
bincountsrhPhi = bincounts(hrhPhi) ./ maximum(bincounts(hrhPhi))
bincountsbbPhi = bincounts(hbbPhi) ./ maximum(bincounts(hbbPhi))
bincountsnu0Phi = bincounts(hnu0Phi) ./ maximum(bincounts(hnu0Phi))


# bincountsm1Phi[end] = 0
# bincountsm2Phi[end] = 0


with_theme(theme_latexfonts()) do
    f = Figure(
        size = (1700, 600), 
        fontsize = 42,  
        figure_padding = 45,
    )
    a = Axis(
        f[1,1], 
        xlabel = L"angular distribution $(^{\circ})$ ",
        ylabel = L"rate (a.u.) $$",
        title = L"theoretical spectra of double beta decay processes $$",
        limits = (0,180.0,0,1.05),
        xticks = (0:90:180, [L"0", L"90", L"180"]),
        yticks = (0:1, [L"0", L"1"]),
    )

    p1 = lines!(a, bcPhi, bincountsnu0Phi, label = L" $0\nu\beta\beta$", linewidth = 8, color = :black)
    p2 = lines!(a, bcPhi, bincountsbbPhi, label = L" $2\nu\beta\beta$", linewidth = 8, color = colors[2], linestyle = (:dot, :dense))
    p3 = lines!(a, bcPhi, bincountsm1Phi, label = L" $0\nu\beta\beta\chi^0, n=1$", linewidth = 8, color = colors[4], linestyle = (:dash, :dense))
    p4 = lines!(a, bcPhi, bincountsm2Phi, label = L" $0\nu\beta\beta\chi^0\chi^0, n=3$", linewidth = 8, color = colors[1], linestyle = (:dashdot, :dense))
    p5 = lines!(a, bc_RH, ys_RH, label = L" $0\nu\beta\beta$ RH", linewidth = 8, color = colors[5], linestyle = (:dashdotdot, :dense))
    # lines!(a, bc, bincountsm2, label = L" $0\nu\beta\beta\chi^2\chi^2$, n=3", linewidth = 8, color = colors[4])
    
    Legend(
        f[1,2], 
        [p1,p2, p3, p4, p5], 
        [L" $0\nu\beta\beta$ ", L" $2\nu\beta\beta$ ", L" $0\nu\beta\beta\chi^0, n=1$ ", L" $0\nu\beta\beta\chi^0\chi^0, n=3$", L" $0\nu\beta\beta$ RH"], 
        patchsize = (80, 50),
        valign = :top,
        width = 400
    )
    # axislegend(a, patchsize = (40, 50), position = :cb, nbanks = 2)
    hidedecorations!(a, label = false,  ticklabels = false, ticks = false,)
    save("scripts/Medex25/figs/theo/phi.png",f, px_per_unit = 5)
    f
end