using DrWatson
@quickactivate "SNSensitivityEstimate"

push!(LOAD_PATH, srcdir())
using ColorSchemes,SensitivityModule, CairoMakie, UnROOT, LaTeXStrings, Revise, FHist, PrettyTables, DataFramesMeta, Combinatorics, LinearAlgebra
using Revise
Revise.track(SensitivityModule)

# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
include(scriptsdir("Params.jl"))

# Load the data from the ROOT files
binsAngle = range(0, 180, 2)
binsESingle = range(0, 3500, 2)
binsESum = range(100, 3500, 35)

all_processes = load_ndim_processes(datadir("sims/fal5_8perc_Boff_TIT_twoDistinct_edepbcu"), binsAngle, binsESingle, binsESum)

processes = [
    get_process("bb0nu_foil_bulk", all_processes),
    get_process("bb_foil_bulk", all_processes),
    get_process("Bi214_foil_bulk", all_processes),
    get_process("Bi214_wire_surface", all_processes),
    get_process("Tl208_foil_bulk", all_processes),
    get_process("K40_foil_bulk", all_processes),
    get_process("Pa234m_foil_bulk", all_processes),
]

set_signal!(processes[2], false)

α = 1.64485362695147
sensitivities = get_sensitivities(SNparams, α, processes; approximate="formula")
maximum(sensitivities)