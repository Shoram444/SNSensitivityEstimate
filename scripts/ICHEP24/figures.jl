using DrWatson

@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, MathTeXEngine, CSV, CairoMakie, UnROOT, ColorSchemes,DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist, Distributions
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

dark_orange = "#EF600A"
light_orange = "#F3DB29"
light_blue = "#2DC3EF"
dark_blue = "#0D3B7C"

### Standard
inPath = joinpath("/home/shoram/Work/PhD_Thesis/Data/Se82/1-Gfactors/spectrumG0.dat")
df = CSV.File(inPath, delim = "    ", header = [" ", "E1", "E2", "dGdE"]) |> DataFrame

df.dGdE ./= maximum(df.dGdE)
gdf = groupby(df, :E1)

# make a projection onto E2 axis
single_spectrum = [maximum(g.dGdE) for g in gdf] 
with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 36,  figure_padding = 36)
    a = Axis(f[1,1], xlabel = L"Single electron energy [keV] $$", ylabel = L"Rate [a.u.] $$" , title = L"Theoretical single-electron energy spectrum $$", )
     
    p = lines!(a, single_spectrum, label = L"standard $2\nu\beta\beta$", color = dark_orange, linewidth = 5)
    
    axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4))


    xlims!(a, 0, 3000)
    ylims!(a, 0, 1.1)

    save(scriptsdir("ICHEP24", "Figs","standard_single.png"), f)
    save(scriptsdir("ICHEP24", "Figs","standard_single.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","standard_single.pdf"), f)
    f
end

df.sumE = round.(df.E1 + df.E2, sigdigits = 4)
gdf2 = groupby(df, :sumE)

sum_spectrum = [sum(g.dGdE) for g in gdf2] 
sum_spectrum ./= maximum(sum_spectrum)

with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 36,  figure_padding = 36)
    a = Axis(f[1,1], xlabel = L"Summed 2-electron energy [keV] $$", ylabel = L"Rate [a.u.] $$" , title = L"Summed 2-electron energy spectrum for $2\nu\beta\beta$", )
     
    p = lines!(a, sum_spectrum, color = dark_orange, linewidth = 5, label= L"standard $2\nu\beta\beta$")

    xlims!(a, 0, 3000)
    ylims!(a, 0, 1.1)
    f
end


using MPGenbb
thetas = [rad2deg(sample_theta_dif(-0.88)) for _=1:1e7]
h = Hist1D(thetas; binedges=0:5:180) 
ys = bincounts(h) ./ maximum(bincounts(h))

with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 36,  figure_padding = 36)
    a = Axis(f[1,1], xlabel = L"Angle between electrons $[^{\circ}]$", ylabel = L"Rate [a.u.] $$" , title = L"Theoretical angular distribution $$", )
     
    p = lines!(a, range(0,180, length(ys)),ys, color = dark_orange, linewidth = 5, label = L"standard $2\nu\beta\beta$")

    xlims!(a, 0, 180)
    ylims!(a, 0, 1.1)
    axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4), position=:cb)
    save(scriptsdir("ICHEP24", "Figs","standard_phi.png"), f)
    save(scriptsdir("ICHEP24", "Figs","standard_phi.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","standard_phi.pdf"), f)
    f
end


### RH
inPath_RH = joinpath("/home/shoram/Work/PhD_Thesis/Data/Se82/Right_Handed/EnDistroRightDataSe10keV.dat")
df_RH = CSV.File(inPath_RH, delim = "\t", header = ["E1", "E2", "dGdE"]) |> DataFrame

df_RH.E1 .*= 1000
df_RH.E2 .*= 1000
df_RH.dGdE ./= maximum(df_RH.dGdE)


gdf_RH = groupby(df_RH, :E1)

# make a projection onto E2 axis
single_spectrum_RH = [maximum(g.dGdE) for g in gdf_RH] 
with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 36,  figure_padding = 36)
    a = Axis(f[1,1], xlabel = L"Single electron energy [keV] $$", ylabel = L"Rate [a.u.] $$" , title = L"Theoretical single-electron energy spectrum $$", )
     
    p = lines!(a, range(-15,3000, length(single_spectrum)),single_spectrum, label = L"standard $2\nu\beta\beta$", color = dark_orange, linewidth = 5)
    lines!(a, range(-15,3000, length(single_spectrum_RH)),single_spectrum_RH, label = L"$2\nu\beta\beta$ with $\bar{\nu}_R$", color = light_blue, linewidth = 5)
    
    axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4))


    xlims!(a, 0, 3000)
    ylims!(a, 0, 1.1)

    save(scriptsdir("ICHEP24", "Figs","RH_single.png"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_single.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_single.pdf"), f)
    f
end

df_RH.sumE = round.(df_RH.E1 + df_RH.E2, sigdigits = 4)
gdf2_RH = groupby(df_RH, :sumE)

sum_spectrum_RH = [sum(g.dGdE) for g in gdf2_RH] 
sum_spectrum_RH ./= maximum(sum_spectrum_RH)

with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 36,  figure_padding = 36)
    a = Axis(f[1,1], xlabel = L"Summed 2-electron energy [keV] $$", ylabel = L"Rate [a.u.] $$" , title = L"Summed 2-electron energy spectrum for $2\nu\beta\beta$", )
     
    p = lines!(a,range(0,3000, length(sum_spectrum)) ,sum_spectrum, color = dark_orange, linewidth = 5, label= L"standard $2\nu\beta\beta$")
    lines!(a,range(0,3000, length(sum_spectrum_RH)) ,sum_spectrum_RH, color = light_blue, linewidth = 5, label= L"$2\nu\beta\beta$ with $\bar{\nu}_R$")

    xlims!(a, 0, 3000)
    ylims!(a, 0, 1.1)
    axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4))
    
    save(scriptsdir("ICHEP24", "Figs","RH_sum.png"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_sum.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_sum.pdf"), f)
    f
end


using MPGenbb
thetas_RH = [rad2deg(sample_theta_dif(0.37)) for _=1:1e7]
h_RH = Hist1D(thetas_RH; binedges=0:5:180) 
ys_RH = bincounts(h_RH) ./ maximum(bincounts(h_RH))

with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 36,  figure_padding = 36)
    a = Axis(f[1,1], xlabel = L"Angle between electrons $[^{\circ}]$", ylabel = L"Rate [a.u.] $$" , title = L"Theoretical angular distribution $$", )
     
    p = lines!(a, range(0,180, length(ys)), ys, color = dark_orange, linewidth = 5, label = L"standard $2\nu\beta\beta$")
    p = lines!(a, range(0,180, length(ys_RH)), ys_RH, color = light_blue, linewidth = 5, label = L"$2\nu\beta\beta$ with $\bar{\nu}_R")

    xlims!(a, 0, 180)
    ylims!(a, 0, 1.1)
    axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4), position=:cb)
    save(scriptsdir("ICHEP24", "Figs","RH_phi.png"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_phi.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_phi.pdf"), f)
    f
end


## Together standard RH

with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 22,  figure_padding = 22)
    Label(f[1, 1:2], "Theoretical predictions", tellwidth = false)
    a1 = Axis(f[3,1],title= L"Energy spectrum $$", xlabel = L"Single electron energy [keV] $$", ylabel = L"Rate [a.u.] $$" )
    a2 = Axis(f[3,2],title= L"Angular distribution $$", xlabel = L"Angle between electrons $[^{\circ}]$", ylabel = L"Rate [a.u.] $$", yaxisposition = :right)
     
    p1 = lines!(a1,range(0,3000, length(sum_spectrum)) ,sum_spectrum, color = dark_orange, linewidth = 3.5, label= L"standard $2\nu\beta\beta$")
    p2 = lines!(a1,range(0,3000, length(sum_spectrum_RH)) ,sum_spectrum_RH, color = light_blue, linewidth = 3.5, label= L"$2\nu\beta\beta$ with $\bar{\nu}_R$")

    p3 = lines!(a2, range(0,180, length(ys)), ys, color = dark_orange, linewidth = 3.5, label = L"standard $2\nu\beta\beta$")
    p4 = lines!(a2, range(0,180, length(ys_RH)), ys_RH, color = light_blue, linewidth = 3.5, label = L"$2\nu\beta\beta$ with $\bar{\nu}_R")

    xlims!(a1, 0, 3000)
    ylims!(a1, 0, 1.1)
    xlims!(a2, 0, 180)
    ylims!(a2, 0, 1.1)

    a2.xticks= ([0,60, 120, 180])
    hideydecorations!(a2, ticks = true, grid = false)
    # axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4), position=:cb)

    Legend(f[2,1:2], [p1, p2], [L"standard $2\nu\beta\beta$",L"$2\nu\beta\beta$ with $\bar{\nu}_R$"], orientation = :horizontal, patchsize = (35, 4))
    colgap!(f.layout, 1, Fixed(38))

    save(scriptsdir("ICHEP24", "Figs","RH_standard_together.png"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_standard_together.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","RH_standard_together.pdf"), f)

    f
end

with_theme(theme_latexfonts()) do 
    f = Figure(size = (750, 500), fontsize = 22,  figure_padding = 22)
    Label(f[1, 1:2], "Theoretical predictions", tellwidth = false)
    a1 = Axis(f[3,1],title= L"Energy spectrum $$", xlabel = L"Single electron energy [keV] $$", ylabel = L"Rate [a.u.] $$" )
    a2 = Axis(f[3,2],title= L"Angular distribution $$", xlabel = L"Angle between electrons $[^{\circ}]$", ylabel = L"Rate [a.u.] $$", yaxisposition = :right)
     
    p1 = lines!(a1,range(0,3000, length(sum_spectrum)) ,sum_spectrum, color = dark_orange, linewidth = 3.5, label= L"standard $2\nu\beta\beta$")
    # p2 = lines!(a1,range(0,3000, length(sum_spectrum_RH)) ,sum_spectrum_RH, color = light_blue, linewidth = 3.5, label= L"$2\nu\beta\beta$ with $\bar{\nu}_R$")

    p3 = lines!(a2, range(0,180, length(ys)), ys, color = dark_orange, linewidth = 3.5, label = L"standard $2\nu\beta\beta$")
    # p4 = lines!(a2, range(0,180, length(ys_RH)), ys_RH, color = light_blue, linewidth = 3.5, label = L"$2\nu\beta\beta$ with $\bar{\nu}_R")

    xlims!(a1, 0, 3000)
    ylims!(a1, 0, 1.1)
    xlims!(a2, 0, 180)
    ylims!(a2, 0, 1.1)

    a2.xticks= ([0,60, 120, 180])
    hideydecorations!(a2, ticks = true, grid = false)
    # axislegend(a, margin = (15, 15, 15, 15), tellwidth = true,patchsize = (35, 4), position=:cb)

    Legend(f[2,1:2], [p1], [L"standard $2\nu\beta\beta$"], orientation = :horizontal, patchsize = (35, 4))
    colgap!(f.layout, 1, Fixed(38))

    save(scriptsdir("ICHEP24", "Figs","standard_together.png"), f)
    save(scriptsdir("ICHEP24", "Figs","standard_together.svg"), f)
    save(scriptsdir("ICHEP24", "Figs","standard_together.pdf"), f)

    f
end

