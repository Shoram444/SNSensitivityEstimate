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