using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie, DataFramesMeta, CSV, Random

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(srcdir("params/Params.jl"))

using UnROOT
f = "data/mva/fal5_TKrec_J40/bb_foil_bulk.root"
vars = [
    "phi", 
    "sumE", 
]
ff = ROOTFile(f)
data = LazyTree(f, "tree", vars) 

bins = (
    phi = (0,180),
    sumE = (300, 3500),
)
roi = (
    phi = (100, 180),
    sumE = (2700, 3100)
)

for (i,n) in enumerate(vars)
    data = filter(x -> getproperty(x, Symbol(n)) > roi[i][1] && getproperty(x, Symbol(n)) < roi[i][2], data)
end

data

using Revise
Revise.track(SensitivityModule)

processes = load_ndim_processes("fal5_TKrec_J40", bins, vars; roi = roi)

signal = get_process("bb0nu_foil_bulk", processes) |> first

hist(getproperty.(signal.data, :sumE))


analysisDict = Dict(
    :Bfield => "Boff", # magnetic field on/off
    :Eres => "8perc", # FWHM of the OMs (sorry for the naming...)
    :mode => "sumE", 
    :trackAlgo => "TKrec",
    :signal => "bb0nu",
    :neutron_config => "current_shielding"
)

# files_directory = "fal5_$(analysisDict[:Eres])_$(analysisDict[:Bfield])_$(analysisDict[:trackAlgo])_twoDistinct_edep_bcu"
files_directory = "mva/fal5_TKrec_J40"

# Load all the processes in the directory. Function `load_processes` takes two arguments:
#  1. dir::String -> the name of the directory where the root files are stored
#  2. mode::String -> the "mode" means which, which dimension we want to investigate, three options (for now) are "sumE", "singleE", "phi"
data_processes = load_data_processes(
    files_directory, 
    "singleE", #analysisDict[:mode],
    fwhm = 0.0,
    roi = roi,
)

sign = get_process("bb_foil_bulk", data_processes) |> first

sign.dataVector

for (field, (minVal, maxVal)) in pairs(roi)
    @show field
    @show minVal, maxVal
end

keys(roi) .|> string |> collect