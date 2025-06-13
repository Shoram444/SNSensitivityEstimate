
###############################################################################
###############################################################################
#### Single energy spectra
###############################################################################
###############################################################################



### Bi214
singleEBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_bulk], 
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleEBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_surface], 
    :bins => binningDict[:SingleE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)


singleEBi214_PMT_glass_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_glass_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

singleEBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

singleEBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

singleEBi214_field_wires_Params = Dict(
    :isotopeName => "Bi214_field_wires", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_field_wires], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_field_wires],
    :bins => binningDict[:SingleE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

singleEBi214_hall_surface_Params = Dict(
    :isotopeName => "Bi214_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_hall_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
)

### Bi210
singleEBi210_foil_bulk_Params = Dict(
    :isotopeName => "Bi210_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_foil_bulk], 
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleEBi210_foil_surface_Params = Dict(
    :isotopeName => "Bi210_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_foil_surface], 
    :bins => binningDict[:SingleE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

singleEBi210_wire_bulk_Params = Dict(
    :isotopeName => "Bi210_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_wire_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

singleEBi210_wire_surface_Params = Dict(
    :isotopeName => "Bi210_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_wire_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

### Tl208

singleETl208_foil_bulk_Params = Dict(
    :isotopeName => "Tl208_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleETl208_foil_surface_Params = Dict(
    :isotopeName => "Tl208_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_foil_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_surface",
    :amount => SNparams["gasVolume"]
)

singleETl208_hall_surface_Params = Dict(
    :isotopeName => "Tl208_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_hall_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
)


singleETl208_PMT_glass_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_glass_bulk], 
    :bins => binningDict[:SingleE],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

### Pa234m

singleEPa234m_foil_bulk_Params = Dict(
    :isotopeName => "Pa234m_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Pa234m_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### K40

singleEK40_foil_bulk_Params = Dict(
    :isotopeName => "K40_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleEK40_hall_surface_Params = Dict(
    :isotopeName => "K40_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_hall_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
)

singleEK40_PMT_glass_bulk_Params = Dict(
    :isotopeName => "K40_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_PMT_glass_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "PMT_glass_bulk",
    :amount => SNparams["PMTGlassMass"]
)

singleEgamma_experimental_surface_Params = Dict(
    :isotopeName => "gamma_experimental_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:gamma_experimental_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:gamma_experimental_surface],
    :bins => binningDict[:SingleE],
    :vertexPosition => "experimental_surface",
    :amount => 1.0
)

### bb_foil_bulk

singleEbbParams = Dict(
    :isotopeName => "bb_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### Xi31
singleEXiParams = Dict(
    :isotopeName => "Xi037_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Xi037_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### RH037
singleERH037Params = Dict(
    :isotopeName => "RH037_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH037_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH037_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### RH020
singleERH020Params = Dict(
    :isotopeName => "RH020_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH020_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH020_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH025
singleERH025Params = Dict(
    :isotopeName => "RH025_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH025_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH025_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH030  
singleERH030Params = Dict(
    :isotopeName => "RH030_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH030_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH030_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH035
singleERH035Params = Dict(
    :isotopeName => "RH035_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH035_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH035_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH040
singleERH040Params = Dict(
    :isotopeName => "RH040_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH040_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH040_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH045
singleERH045Params = Dict(
    :isotopeName => "RH045_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH045_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH045_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
### RH050
singleERH050Params = Dict(
    :isotopeName => "RH050_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH050_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH050_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)


### bb0nu_foil_bulk

singleEbb0nuParams = Dict(
    :isotopeName => "bb0nu_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nu_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nu_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### bb0nuM1_foil_bulk

singleEbb0nuM1Params = Dict(
    :isotopeName => "bb0nuM1_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuM1_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuM1_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### bb0nuM2_foil_bulk

singleEbb0nuM2Params = Dict(
    :isotopeName => "bb0nuM2_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuM2_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuM2_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### bb0nuRHl_foil_bulk

singleEbb0nuRHlParams = Dict(
    :isotopeName => "bb0nuRHl_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuRHl_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuRHl_foil_bulk],
    :bins => binningDict[:SingleE],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)


### neutron_external

singleNeutronParams = Dict(
    :isotopeName => "neutron_external",
    :signal => :false,
    :bins => binningDict[:SingleE],
)

singleEParams = Dict(
    :Bi214_foil_bulk => singleEBi214_foil_bulk_Params,
    :Bi214_foil_surface => singleEBi214_foil_surface_Params,
    :Bi214_PMT_glass_bulk => singleEBi214_PMT_glass_bulk_Params,
    :Bi214_wire_bulk => singleEBi214_wire_bulk_Params,
    :Bi214_wire_surface => singleEBi214_wire_surface_Params,
    :Bi214_field_wires => singleEBi214_field_wires_Params,
    :Bi214_hall_surface => singleEBi214_hall_surface_Params,
    :Bi210_foil_bulk => singleEBi210_foil_bulk_Params,
    :Bi210_foil_surface => singleEBi210_foil_surface_Params,
    :Bi210_wire_bulk => singleEBi210_wire_bulk_Params,
    :Bi210_wire_surface => singleEBi210_wire_surface_Params,
    :Tl208_foil_bulk => singleETl208_foil_bulk_Params,
    :Tl208_foil_surface => singleETl208_foil_surface_Params,
    :Tl208_hall_surface => singleETl208_hall_surface_Params,
    :Tl208_PMT_glass_bulk => singleETl208_PMT_glass_bulk_Params,
    :Pa234m_foil_bulk => singleEPa234m_foil_bulk_Params,
    :K40_foil_bulk => singleEK40_foil_bulk_Params,
    :K40_PMT_glass_bulk => singleEK40_PMT_glass_bulk_Params,
    :K40_hall_surface => singleEK40_hall_surface_Params,
    :bb_foil_bulk => singleEbbParams,
    :Xi037_foil_bulk => singleEXiParams,
    :RH037_foil_bulk => singleERH037Params,
    :RH020_foil_bulk => singleERH020Params,
    :RH025_foil_bulk => singleERH025Params,
    :RH030_foil_bulk => singleERH030Params,
    :RH035_foil_bulk => singleERH035Params,
    :RH040_foil_bulk => singleERH040Params,
    :RH045_foil_bulk => singleERH045Params,
    :RH050_foil_bulk => singleERH050Params,
    :bb0nu_foil_bulk => singleEbb0nuParams,
    :bb0nuM1_foil_bulk => singleEbb0nuM1Params,
    :bb0nuM2_foil_bulk => singleEbb0nuM2Params,
    :bb0nuRHl_foil_bulk => singleEbb0nuRHlParams,
    :neutron_external => singleNeutronParams,
    :gamma_experimental_surface => singleEgamma_experimental_surface_Params
)
