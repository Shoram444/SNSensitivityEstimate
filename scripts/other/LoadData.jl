using DrWatson
@quickactivate "SensitivityEstimate"

using Revise
push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"])

function load_files(dir::String)
    filesDict = Dict()
    full_dir = datadir("sims", dir)

    println("Loading files from: $full_dir ...")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        df = ROOTFile(joinpath(full_dir, file)) |> ffrf
        fileName = split(file, ".")[1]  |> split |> first 

        filesDict[fileName] = df
    end
    println("Loaded $nFiles files.")

    return filesDict
end

function load_processes(dir::String, mode::String)
    filesDict = Dict()
    full_dir = datadir("sims", dir)
    processes = []

    println("Loading files from: $full_dir ...")
    println("mode: $mode ")

    nFiles = 0
    for file in readdir(full_dir)
        nFiles += 1
        df = ROOTFile(joinpath(full_dir, file)) |> ffrf
        fileName = split(file, ".")[1]  |> split |> first 
        
        if mode == "sumE"
            push!(
                processes,
                Process(
                    df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
                    sumEParams[Symbol(fileName)]
                )
            )
        elseif mode == "singleE"
            push!(
                processes,
                Process(
                    vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
                    singleEParams[Symbol(fileName)]
                )
            )
        elseif mode == "phi"
            push!(
                processes,
                Process(
                    df.phi,
                    phiParams[Symbol(fileName)]
                )
            )
        else
            @error "mode must be one of: sumE, singleE, phi. Chosen mode: $mode!"
        end
    end
    println("Loaded $nFiles files.")

    return processes
end


####################################################
###### SDBDRC FILES
####################################################
# SDBDRCFiles = Dict()
# SDBDRCDirectory = datadir("sims", "SDBDRC/")

# for file in readdir(SDBDRCDirectory)
#     df = ROOTFile(joinpath(SDBDRCDirectory, file)) |> ffrf
#     fileName = join(split(file, "_")[1:end-1], "_")

#     SDBDRCFiles[fileName] = df
# end

# ####################################################
# ###### SDBDRC+vertex Cut FILES
# ####################################################
# vertexFiles = Dict()
# vertexDirectory = datadir("sims", "SDBDRC_vertex/")

# for file in readdir(vertexDirectory)
#     df = ROOTFile(joinpath(vertexDirectory, file)) |> ffrf
#     fileName = join(split(file, "_")[1:end-3], "_")

#     vertexFiles[fileName] = df
# end

# ####################################################
# ###### SDBDRC+vertex+Probability Cut FILES
# ####################################################
# probFiles = Dict()
# probProcessesESum = []
# probProcessesESumMeV = []
# probProcessesESingle = []
# probProcessesPhi = []
# probDirectory = datadir("sims", "SDBDRC_vertex_prob/")

# @show probDirectory
# for file in readdir(probDirectory)
#     df = ROOTFile(joinpath(probDirectory, file)) |> ffrf
#     fileName = join(split(file, "_")[1:end-3], "_")

#     probFiles[fileName] = df
#     push!(
#         probProcessesESum,
#         Process(
#             df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
#             sumEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesESumMeV,
#         Process(
#             (df.reconstructedEnergy1 .+ df.reconstructedEnergy2) ./ 1000.0,
#             sumEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesESingle,
#         Process(
#             vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
#             singleEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesPhi,
#         Process(
#             df.phi,
#             phiParams[Symbol(fileName)]
#         )
#     )
# end

# ####################################################
# ###### 11% resolution
# ###### SDBDRC+vertex+Probability Cut FILES
# ####################################################
# probFiles11 = Dict()
# probProcessesESum11 = []
# probProcessesESingle11 = []
# probProcessesPhi11 = []
# probDirectory11 = datadir("sims", "prob_11percent_res/")

# @show probDirectory11
# for file in readdir(probDirectory11)
#     df = ROOTFile(joinpath(probDirectory11, file)) |> ffrf
#     fname = split(file, "_")
#     fileName = join(vcat(fname[1], fname[3:end-1]), "_")

#     probFiles11[fileName] = df
#     push!(
#         probProcessesESum11,
#         Process(
#             df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
#             sumEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesESingle11,
#         Process(
#             vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
#             singleEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesPhi11,
#         Process(
#             df.phi,
#             phiParams[Symbol(fileName)]
#         )
#     )
# end

# ####################################################
# ###### Fal 5 8%
# ###### SDBDRC+vertex+Probability Cut FILES
# ####################################################
# probFiles11 = Dict()
# probProcessesESum11 = []
# probProcessesESingle11 = []
# probProcessesPhi11 = []
# probDirectory11 = datadir("sims", "prob_11percent_res/")

# @show probDirectory11
# for file in readdir(probDirectory11)
#     df = ROOTFile(joinpath(probDirectory11, file)) |> ffrf
#     fname = split(file, "_")
#     fileName = join(vcat(fname[1], fname[3:end-1]), "_")

#     probFiles11[fileName] = df
#     push!(
#         probProcessesESum11,
#         Process(
#             df.reconstructedEnergy1 .+ df.reconstructedEnergy2,
#             sumEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesESingle11,
#         Process(
#             vcat(df.reconstructedEnergy1, df.reconstructedEnergy2),
#             singleEParams[Symbol(fileName)]
#         )
#     )

#     push!(
#         probProcessesPhi11,
#         Process(
#             df.phi,
#             phiParams[Symbol(fileName)]
#         )
#     )
# end

