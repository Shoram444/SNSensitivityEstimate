using DrWatson
@quickactivate "SensitivityEstimate"

push!(LOAD_PATH, srcdir())
using SensitivityModule, StatsPlots, UnROOT, DataFramesMeta, LaTeXStrings, Revise, StatsBase, FHist
Revise.track(SensitivityModule)

include(scriptsdir("Params.jl"))

###################################################
###################################################
#### LOADING FILES
###################################################
###################################################

# hack for quick fill from root
ffrf(file) = fill_from_root_file(file, "tree", ["phi", "reconstructedEnergy1", "reconstructedEnergy2"] ) 

#Bi214
Bi214_foil_bulk = ROOTFile("data/sims/Bi214_foil_bulk_1e8E.root") |> ffrf
Bi214_foil_surface = ROOTFile("data/sims/Bi214_foil_surface_1e7E.root") |> ffrf
Bi214_PMT_bulk = ROOTFile("data/sims/Bi214_PMT_bulk_1e7E.root") |> ffrf
Bi214_wire_bulk = ROOTFile("data/sims/Bi214_wire_bulk_1e7E.root") |> ffrf
Bi214_wire_surface = ROOTFile("data/sims/Bi214_wire_surface_1e7E.root") |> ffrf

#Tl208
Tl208_foil_bulk = ROOTFile("data/sims/Tl208_foil_bulk_1e8E.root") |> ffrf
Tl208_PMT_bulk = ROOTFile("data/sims/Tl208_PMT_bulk_1e8E.root") |> ffrf

#Pa234m
Pa234m_foil_bulk = ROOTFile("data/sims/Pa234m_foil_bulk_1e8E.root") |> ffrf

#K40
K40_foil_bulk = ROOTFile("data/sims/K40_foil_bulk_1e8E.root") |> ffrf
K40_PMT_bulk = ROOTFile("data/sims/K40_PMT_bulk_1e8E.root") |> ffrf

#Xi31
Xi31_foil_bulk = ROOTFile("data/sims/Xi037_foil_bulk_1e8E.root") |> ffrf


###################################################
###################################################
#### INITIATE PROCESS ESUM
###################################################
###################################################

Bi214_foil_bulk_SumE    = Process(Bi214_foil_bulk.reconstructedEnergy1 .+ Bi214_foil_bulk.reconstructedEnergy2, sumEBi214_foil_bulk_Params)
Bi214_foil_surface_SumE = Process(Bi214_foil_surface.reconstructedEnergy1 .+ Bi214_foil_surface.reconstructedEnergy2, sumEBi214_foil_surface_Params)  
Bi214_PMT_bulk_SumE     = Process(Bi214_PMT_bulk.reconstructedEnergy1 .+ Bi214_PMT_bulk.reconstructedEnergy2, sumEBi214_PMT_bulk_Params) 
Bi214_wire_bulk_SumE    = Process(Bi214_wire_bulk.reconstructedEnergy1 .+ Bi214_wire_bulk.reconstructedEnergy2, sumEBi214_wire_bulk_Params) 
Bi214_wire_surface_SumE = Process(Bi214_wire_surface.reconstructedEnergy1 .+ Bi214_wire_surface.reconstructedEnergy2, sumEBi214_wire_surface_Params) 

Tl208_foil_bulk_SumE    = Process(Tl208_foil_bulk.reconstructedEnergy1 .+ Tl208_foil_bulk.reconstructedEnergy2, sumETl208_foil_bulk_Params)
Tl208_PMT_bulk_SumE     = Process(Tl208_PMT_bulk.reconstructedEnergy1 .+ Tl208_PMT_bulk.reconstructedEnergy2, sumETl208_PMT_bulk_Params)

