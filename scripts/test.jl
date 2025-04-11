import Pkg
Pkg.activate("/pbs/home/m/mpetro/sps_mpetro/Projects/PhD/SNSensitivityEstimate")
Pkg.instantiate()
# println("start")
using DrWatson
@quickactivate "SNSensitivityEstimate"
println("2")

push!(LOAD_PATH, srcdir())
println("3")

using SensitivityModule
println("4")

using Random, LinearAlgebra, Statistics, Distributions, Plots, BAT, BinnedModels, StatsBase, DensityInterface, IntervalSets, FHist, SpecialFunctions, ValueShapes
global t0 = time()

analysisDict = Dict(
    :signal => "bb0nu_foil_bulk",
    :roi => (300, 3400),
    :bw => 100,
    :mode => "sumE"
)

println("great success!")
