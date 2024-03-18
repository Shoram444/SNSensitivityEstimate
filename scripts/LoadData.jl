using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])

####################################################
###### SDBDRC FILES
####################################################
SDBDRCFiles = Dict()
SDBDRCDirectory = datadir("sims", "SDBDRC/")

for file in readdir(SDBDRCDirectory)
    df = ROOTFile(joinpath(SDBDRCDirectory, file)) |> ffrf
    fileName = join(split(file, "_")[1:end-1], "_")

    SDBDRCFiles[fileName] = df
end

####################################################
###### SDBDRC+vertex Cut FILES
####################################################
vertexFiles = Dict()
vertexDirectory = datadir("sims", "SDBDRC_vertex/")

for file in readdir(vertexDirectory)
    df = ROOTFile(joinpath(vertexDirectory, file)) |> ffrf
    fileName = join(split(file, "_")[1:end-3], "_")

    vertexFiles[fileName] = df
end

####################################################
###### SDBDRC+vertex+Probability Cut FILES
####################################################
probFiles = Dict()
probProcessesESum = []
probProcessesESingle = []
probProcessesPhi = []
probDirectory = datadir("sims", "SDBDRC_vertex_prob/")

for file in readdir(probDirectory)
    df = ROOTFile(joinpath(probDirectory, file)) |> ffrf
    fileName = join(split(file, "_")[1:end-3], "_")

    probFiles[fileName] = df
    push!(
        probProcessesESum,
        Process(
            df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
            sumEParams[Symbol(fileName)]
        )
    )

    push!(
        probProcessesESingle,
        Process(
            vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
            singleEParams[Symbol(fileName)]
        )
    )

    push!(
        probProcessesPhi,
        Process(
            df.phi,
            phiParams[Symbol(fileName)]
        )
    )
end

