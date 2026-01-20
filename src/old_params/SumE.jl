#### Sum energy spectra


### Bi214
sumEBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_bulk], 
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumEBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_surface], 
    :bins => binningDict[:SumE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

sumEBi214_PMT_glass_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_glass_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

sumEBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

sumEBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

sumEBi214_field_wires_Params = Dict(
    :isotopeName => "Bi214_field_wires", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_field_wires], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_field_wires],
    :bins => binningDict[:SumE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

sumEBi214_hall_surface_Params = Dict(
    :isotopeName => "Bi214_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_hall_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
)

### Bi210
sumEBi210_foil_bulk_Params = Dict(
    :isotopeName => "Bi210_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_foil_bulk], 
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumEBi210_foil_surface_Params = Dict(
    :isotopeName => "Bi210_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_foil_surface], 
    :bins => binningDict[:SumE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"],
)


sumEBi210_wire_bulk_Params = Dict(
    :isotopeName => "Bi210_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_wire_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

sumEBi210_wire_surface_Params = Dict(
    :isotopeName => "Bi210_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_wire_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

### Tl208

sumETl208_foil_bulk_Params = Dict(
    :isotopeName => "Tl208_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumETl208_foil_surface_Params = Dict(
    :isotopeName => "Tl208_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_foil_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_surface",
    :amount => SNparams["gasVolume"]
)

sumETl208_hall_surface_Params = Dict(
    :isotopeName => "Tl208_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_hall_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
)

sumETl208_PMT_glass_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_glass_bulk], 
    :bins => binningDict[:SumE],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

### Pa234m

sumEPa234m_foil_bulk_Params = Dict(
    :isotopeName => "Pa234m_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Pa234m_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### K40

sumEK40_foil_bulk_Params = Dict(
    :isotopeName => "K40_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumEK40_hall_surface_Params = Dict(
    :isotopeName => "K40_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_hall_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
)

sumEK40_PMT_glass_bulk_Params = Dict(
    :isotopeName => "K40_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_PMT_glass_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "PMT_glass_bulk",
    :amount => SNparams["PMTGlassMass"]
)

sumEgamma_experimental_surface_Params = Dict(
    :isotopeName => "gamma_experimental_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:gamma_experimental_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:gamma_experimental_surface],
    :bins => binningDict[:SumE],
    :vertexPosition => "experimental_surface",
    :amount => 1.0
)

sumEgamma_hall_bulk_Params = Dict(
    :isotopeName => "gamma_hall_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:gamma_hall_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:gamma_hall_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "hall_bulk",
    :amount => 1.0
)

### 2nubb

sumEbbParams = Dict(
    :isotopeName => "bb_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### Xi31
sumEXiParams = Dict(
    :isotopeName => "Xi037_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Xi037_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### RH037
sumERH037Params = Dict(
    :isotopeName => "RH037_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH037_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH037_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### RH020
sumERH020Params = Dict(
    :isotopeName => "RH020_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH020_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH020_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH025
sumERH025Params = Dict(
    :isotopeName => "RH025_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH025_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH025_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH030
sumERH030Params = Dict(
    :isotopeName => "RH030_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH030_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH030_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH035
sumERH035Params = Dict(
    :isotopeName => "RH035_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH035_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH035_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH040
sumERH040Params = Dict(
    :isotopeName => "RH040_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH040_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH040_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH045
sumERH045Params = Dict(
    :isotopeName => "RH045_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH045_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH045_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH050
sumERH050Params = Dict(
    :isotopeName => "RH050_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH050_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH050_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubb
sumEbb0nuParams = Dict(
    :isotopeName => "bb0nu_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nu_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nu_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubbM1
sumEbb0nuM1Params = Dict(
    :isotopeName => "bb0nuM1_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuM1_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuM1_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubbM2
sumEbb0nuM2Params = Dict(
    :isotopeName => "bb0nuM2_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuM2_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuM2_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubbScalar0
sumEbb0nuScalar0Params = Dict(
    :isotopeName => "bb0nuScalar0_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuScalar0_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuScalar0_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubbScalar1
sumEbb0nuScalar1Params = Dict(
    :isotopeName => "bb0nuScalar1_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuScalar1_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuScalar1_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubbScalar2
sumEbb0nuScalar2Params = Dict(
    :isotopeName => "bb0nuScalar2_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuScalar2_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuScalar2_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### 0nubbRHl
sumEbb0nuRHlParams = Dict(
    :isotopeName => "bb0nuRHl_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuRHl_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuRHl_foil_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### neutron
sumNeutronParams = Dict(
    :isotopeName => "neutron_external", 
    :signal => :false, 
    :bins => binningDict[:SumE],
)

### K40 8inch calo bulk
sumEK40_calo_8inch_scin_bulk_Params = Dict(
    :isotopeName => "K40_calo_8inch_scin_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_calo_8inch_scin_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_calo_8inch_scin_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "calo_8inch_scin_bulk",
    :amount => SNparams["caloMass_8inch"]
)

### Bi214 8inch calo bulk
sumEBi214_calo_8inch_scin_bulk_Params = Dict(
    :isotopeName => "Bi214_calo_8inch_scin_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_calo_8inch_scin_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_calo_8inch_scin_bulk],
    :bins => binningDict[:SumE],
    :vertexPosition => "calo_8inch_scin_bulk",
    :amount => SNparams["caloMass_8inch"]
)

sumEParams = Dict(
    :Bi214_foil_bulk => sumEBi214_foil_bulk_Params,
    :Bi214_foil_surface => sumEBi214_foil_surface_Params,
    :Bi214_hall_surface => sumEBi214_hall_surface_Params,
    :Bi214_PMT_glass_bulk => sumEBi214_PMT_glass_bulk_Params,
    :Bi214_wire_bulk => sumEBi214_wire_bulk_Params,
    :Bi214_wire_surface => sumEBi214_wire_surface_Params,
    :Bi214_field_wires => sumEBi214_field_wires_Params,
    :Bi210_foil_bulk => sumEBi210_foil_bulk_Params,
    :Bi210_foil_surface => sumEBi210_foil_surface_Params,
    :Bi210_wire_bulk => sumEBi210_wire_bulk_Params,
    :Bi210_wire_surface => sumEBi210_wire_surface_Params,
    :Tl208_foil_bulk => sumETl208_foil_bulk_Params,
    :Tl208_foil_surface => sumETl208_foil_surface_Params,
    :Tl208_hall_surface => sumETl208_hall_surface_Params,
    :Tl208_PMT_glass_bulk => sumETl208_PMT_glass_bulk_Params,
    :Pa234m_foil_bulk => sumEPa234m_foil_bulk_Params,
    :K40_foil_bulk => sumEK40_foil_bulk_Params,
    :K40_hall_surface => sumEK40_hall_surface_Params,
    :K40_PMT_glass_bulk => sumEK40_PMT_glass_bulk_Params,
    :bb_foil_bulk => sumEbbParams,
    :Xi037_foil_bulk => sumEXiParams,
    :RH037_foil_bulk => sumERH037Params,
    :RH020_foil_bulk => sumERH020Params,
    :RH025_foil_bulk => sumERH025Params,
    :RH030_foil_bulk => sumERH030Params,
    :RH035_foil_bulk => sumERH035Params,
    :RH040_foil_bulk => sumERH040Params,
    :RH045_foil_bulk => sumERH045Params,
    :RH050_foil_bulk => sumERH050Params,
    :bb0nu_foil_bulk => sumEbb0nuParams,
    :bb0nuM1_foil_bulk => sumEbb0nuM1Params,
    :bb0nuM2_foil_bulk => sumEbb0nuM2Params,
    :bb0nuRHl_foil_bulk => sumEbb0nuRHlParams,
    :neutron_external => sumNeutronParams,
    :gamma_experimental_surface => sumEgamma_experimental_surface_Params,
    :gamma_hall_bulk => sumEgamma_hall_bulk_Params,
    :K40_calo_8inch_scin_bulk => sumEK40_calo_8inch_scin_bulk_Params,
    :Bi214_calo_8inch_scin_bulk => sumEBi214_calo_8inch_scin_bulk_Params,
    :bb0nuScalar0_foil_bulk => sumEbb0nuScalar0Params,
    :bb0nuScalar1_foil_bulk => sumEbb0nuScalar1Params,
    :bb0nuScalar2_foil_bulk => sumEbb0nuScalar2Params
)