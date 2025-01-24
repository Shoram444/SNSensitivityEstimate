# SNSensitivityEstimate

All **performed and completed tasks** are summarized in the [TO DO](https://github.com/Shoram444/SNSensitivityEstimate/issues/14).

This code repository represents a small project within the SuperNEMO experiment. The goal is calcualtion of sensitivity estimates for various modes of $2\nu\beta\beta$ decay. 

The data is available at `/sps/nemo/scratch/mpetro/Projects/SNSensitivityEstimate/data` and should be copied to the `data/` directory. 

The methodology of calculating sensitivities is described in `docdb:4816` with a short presentation given in `docdb:5568` and **more recently** `docdb:5833`. The `plots/` folder in this repository contains all the relevant figures for finding ROI and shows the individual spectra.

**EXAMPLES**
========

Main example of how to use this repository to calculate sensitivity to 0nu is provided in `scripts/0nu/Example_bb0nu_Sensitivity.jl`.

There are example notebooks in the `notebooks/` folder. 

**DATA**
========
## Input data format for this repository:
The input data **must** be provided in form of a `root` file, which contains (at least) the following information **per event**:
- `float reconstructedEnergy1`
- `float reconstructedEnergy2`
- `float phi`

Each filename must match the naming convention: `isotope_vertex_generator.root`, (i.e. `Bi214_foil_bulk.root`) for the process to be compatible with the `LoadData::load_processes(dir::String, mode::String)` function. 
( Optionally you can initiaite your own processes and work with them manually, see documentation in `src/Process.jl` or `? Process` in `REPL`. )

A very important file is the: `scripts/Params.jl`. This contains the dictionaries for all the processes/acitvitites/info about the detector/etc. Look through the file and make sure you follow the conventions there. 

## Most recent simulation performed is in `/sps/nemo/scratch/mpetro/Projects/Phd/SNSensitivityEstimate/data/sims/fal5_12perc_Boff_TIT`:

## Simulated processes
- [x] 2nubb - foil_bulk
- [x] 0nubb - foil_bulk
- [x] Bi214 - foil_bulk
- [x] Tl208 - foil_bulk
- [ ] Pa234m - foil_bulk
- [ ] K40 - foil_bulk
- [x] Bi214 - wire_surface (radon)

## Data-cuts using `SNCuts`
```bash
useEventHasTwoTracks : boolean = true
useEventHasTwoFoilVertices : boolean = true
useEventHasTwoCaloHits : boolean = true
useEventHasTwoDistinctAssociatedCaloHits : boolean = true

useEventHasSumEnergyAbove : boolean = true
minSumEnergy : real = 300.0
useEventHasSumEnergyBelow : boolean = true
maxSumEnergy : real = 3500.0

useEventHasFoilVertexDistanceBelow : boolean = true
maxFoilVertexDistance : real = 50.0  # mm

useEventHasPintAbove : boolean = true  # internal probability ToF cut
minPint : real = 0.04

useEventHasPextBelow : boolean = true   # external probability ToF cut
maxPext : real = 0.01
```

## Simulation setup
- Use `realistic_flat` source foil geometry 
- Use magnetic field off
- Use TKReconstrcut to reconstruct the trajectories (with polylines)
- Use 12% FWHM across main OMs (Or do we have some better solution ready?)
- Use basic SimRC setup


**CODE**
========

The main bulk of the analysis is done in the file `scripts/main.jl`. There is a `julia module` provided within this repository in the `src` folder. There are a number of small scripts in the `scripts` folder which mainly deal with plotting and miscelnaceous. 

**REPRODUCIBILITY**
========
## To get Julia
**Download Juliaup**
```bash
curl -fsSL https://install.julialang.org | sh
``` 
**install version 1.10.0**
```bash
juliaup add 1.10.0
```


## To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data not included in the git repository.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```


This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "SNSensitivityEstimate"
```
which auto-activate the project and enable local path handling from DrWatson.
