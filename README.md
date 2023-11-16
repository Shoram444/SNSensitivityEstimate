# SNSensitivityEstimate

This code repository represents a small project within the SuperNEMO experiment. The goal is calcualtion of sensitivity estimates for various modes of $2\nu\beta\beta$ decay. 

The data is available at `/sps/nemo/scratch/mpetro/Projects/SNSensitivityEstimate/data` and should be copied to the `data/` directory. 

The methodology of calculating sensitivities is described in `docdb:4816` with a short presentation given in `docdb:5568`. The `plots/` folder in this repository contains all the relevant figures for finding ROI and shows the individual spectra.

**Data**
The data used as of 16th November is described by the following:
1. $10^8$ events were generated for each of the following processes: 
   1. Internal Background processes: ${}^{214}Bi$, ${}^{208}Tl$, ${}^{234m}Pa$, ${}^{40}K$ - the activities of these is obtained from table 1 in [paper](https://link.springer.com/article/10.1140/epjc/s10052-018-6295-x)
   2. Signal processes: "standard" $2\nu\beta\beta$ 
2. The following data-cuts were applied to data:
   1. two negatively charged particle tracks reconstructed,
   2. two vertices on the source foil,
   3. sum of reconstructed electron energies within the range: $(0 - 3500) keV$,
   4. two individual OM hits,
   5. two associated OM hits. 
3. The following variables are extracted from the simulation:
   1. Individual electron energies (2 per event)
   2. The angle of escape between the two electrons - obtained by drawing a line between the vertex on foil and on the calorimeter **NOT FROM TRACKING ALGORITHM!!**
4. The so-called *mock* calibration was used (this is roughly 8% at 1MeV) to obtain reconstructed energies

**Code**
The main bulk of the analysis is done in the file `main.jl`. There is a `julia module` provided within this repository in the `src` folder. There are a number of small scripts in the `scripts` folder which mainly deal with plotting and miscelnaceous. 

**To (locally) reproduce this project, do the following**:

0. Download this code base. Notice that raw data are typically not included in the git-history and may need to be downloaded independently.
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