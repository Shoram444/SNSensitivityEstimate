using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
using PrettyTables
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))
include(scriptsdir("LoadData.jl"))

signalProcessName = "bb0nu_foil_bulk"
sigProcessESum = get_process(signalProcessName, probProcessesESum)


bb2nuProcess = get_process("bb_foil_bulk", probProcessesESum)
bb2nuProcess.signal = false
bkgProcessesESum = [
    get_process("Tl208_foil_bulk", probProcessesESum),
    get_process("Bi214_foil_bulk", probProcessesESum),
    get_process("Bi214_foil_surface", probProcessesESum),
    get_process("Bi214_wire_surface", probProcessesESum),
    bb2nuProcess
]


nExpBkg = Float64[]
sigEff = Float64[]
bkgEff = Float64[]
bkgIsotopes = String[]

# for p in bkgProcessesESum
#     push!( bkgIsotopes, p.isotopeName )
# end

minROI = [2500, 2600, 2700, 2700, 2800] .+ 50
maxROI = [3000, 3200, 3100, 3200, 3200] .- 50

ROIs = []

for (min, max) in zip(minROI, maxROI)
    for p in bkgProcessesESum
        push!(bkgIsotopes, p.isotopeName)
        push!(ROIs, (min - 50, max +50))
        push!(nExpBkg, get_bkg_counts_ROI(min, max, p))
        push!(sigEff, lookup(sigProcessESum, min, max))
        push!(bkgEff, lookup(p, min, max))
    end
end

dfbkg = DataFrame(
    ROI=ROIs,
    isotope=bkgIsotopes,
    nb=nExpBkg,
    ε=sigEff .* 100,
    εbkg=bkgEff .* 100
)

gdf = groupby(dfbkg, :ROI)

for g in gdf
    println("ROI: $(g[1,1]) keV; εSignal = $(g[1,end-1]) %\n")

    pretty_table(
        select(g, Not(:ROI, :ε)),
        header=["Bkgprocess", "nb", "effBkg"],
        backend=Val(:markdown)
    )
    println()
end
