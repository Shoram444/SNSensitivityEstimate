# include(srcdir("SensitivityModule.jl")) # to get halfLife_to_activity
import SensitivityModule.halfLife_to_activity
using Measurements
#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23,                         # Avogadro's number in [1/mol]
    "W" => 0.08192,                             # Se82 molar mass in [kg/mol]
    "a" => 0.98,                                # abundance/foil enrichment; check number
    "foilMass" => 6.25,                         # foil mass in [kg]
    "gasVolume" => 15,                          # tracker volume in [m3]
    "PMTGlassMass" => 286,                      # PMT glass mass in [kg]
    "wireBulkMass" => 20,                       # DUMMY VALUE!! mass of the tracker wires in [kg]
    "t" => 2.5 * 365 * 24 * 3600,               # measurement time in [s]
    "tYear" => 2.5,                             # measurement time in [y]
    "SeThalf2nu" => (9.39 ± 0.17)*1e19 * 365 * 24 * 3600,    # 2nu Se82 half life in [s], results from NEMO-3
    "SeThalf0nu" => 1e26 * 365 * 24 * 3600,    # 0nu Se82 half life in [s], results from NEMO-3
) 

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x and docDB 4505
    :Pa234m_foil_bulk => (17.3 ± 0.1) / 1000 ,          # [mBq/kg] converted to [Bq/kg] from NEMO-3
    :Bi214_foil_bulk => 10 / 1_000_000 ,               # [μBq/kg] converted to [Bq/kg] from SN measurements
    :Bi214_foil_surface => 150 / 1_000_000 ,                  # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3) 
    :Bi214_wire_surface => 150 / 1_000_000 ,                  # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3) 
    :Bi214_wire_bulk => 0.00001 ,                       # MOCK VALUE
    :Bi214_PMT_glass_bulk => 140 / 286,                 # [Bq/kg] originally the value is given as 417Bq, I just divide by PMT weight here
    :Tl208_foil_bulk => 2 / 1_000_000 ,                # [μBq/kg] converted to [Bq/kg] from SN measurements
    :Tl208_foil_surface => 150 / 1_000_000 ,                  # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3)
    :Tl208_PMT_glass_bulk => 41.4 / 286 ,               # [Bq/kg] originally the value is given as 41.4Bq, I just divide by PMT weight here
    :Bi210_foil_bulk => 10 / 1_000_000 ,                # MOCK VALUE
    :Bi210_foil_surface => 150 / 1_000_000 ,                  # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3) 
    :Bi210_wire_surface => 150 / 1_000_000 ,                  # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3) 
    :Bi210_wire_bulk => 0.00001 ,                       # MOCK VALUE
    :K40_foil_bulk => (58.7 ± 0.1) / 1000 ,             # [mBq/kg] converted to [Bq/kg] from NEMO-3
    :K40_PMT_glass_bulk => 417 / 286 ,                  # [Bq/kg] originally the value is given as 417Bq, I just divide by PMT weight here
)

SigActivityParams = Dict( 
    :bb_foil_bulk => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # activity calculated from 2nubb Half-life in [Bq/kg]
    :Xi037_foil_bulk => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH037_foil_bulk => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :bb0nu_foil_bulk => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
)

SimulationParams = Dict(
    :Bi214_foil_bulk => 1e8,
    :Bi214_foil_surface => 1e7,
    :Bi214_PMT_glass_bulk => 1e7,
    :Bi214_wire_bulk => 1e7,
    :Bi214_wire_surface => 1e7,
    :Bi210_foil_bulk => 1e6 ,       
    :Bi210_foil_surface => 1e6 ,    
    :Bi210_wire_surface => 1e6 ,    
    :Bi210_wire_bulk => 1e6 ,       
    :Tl208_foil_bulk => 1e8,
    :Tl208_foil_surface => 117e6,
    :Tl208_PMT_glass_bulk => 1e8,
    :Pa234m_foil_bulk => 1e8,
    :K40_foil_bulk => 1e8,
    :K40_PMT_glass_bulk => 1e8,
    :bb_foil_bulk => 1e8,
    :Xi037_foil_bulk => 1e8 ,
    :RH037_foil_bulk => 1e7 ,
    :bb0nu_foil_bulk => 1e7,
)

binningDict = Dict(
    :SumE => 0:100:3500,
    :SingleE => 0:100:3500,
    :Phi => 0:5:180,
)

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

sumEParams = Dict(
    :Bi214_foil_bulk => sumEBi214_foil_bulk_Params,
    :Bi214_foil_surface => sumEBi214_foil_surface_Params,
    :Bi214_PMT_glass_bulk => sumEBi214_PMT_glass_bulk_Params,
    :Bi214_wire_bulk => sumEBi214_wire_bulk_Params,
    :Bi214_wire_surface => sumEBi214_wire_surface_Params,
    :Bi210_foil_bulk => sumEBi210_foil_bulk_Params,
    :Bi210_foil_surface => sumEBi210_foil_surface_Params,
    :Bi210_wire_bulk => sumEBi210_wire_bulk_Params,
    :Bi210_wire_surface => sumEBi210_wire_surface_Params,
    :Tl208_foil_bulk => sumETl208_foil_bulk_Params,
    :Tl208_foil_surface => sumETl208_foil_surface_Params,
    :Tl208_PMT_glass_bulk => sumETl208_PMT_glass_bulk_Params,
    :Pa234m_foil_bulk => sumEPa234m_foil_bulk_Params,
    :K40_foil_bulk => sumEK40_foil_bulk_Params,
    :K40_PMT_glass_bulk => sumEK40_PMT_glass_bulk_Params,
    :bb_foil_bulk => sumEbbParams,
    :Xi037_foil_bulk => sumEXiParams,
    :RH037_foil_bulk => sumERH037Params,
    :bb0nu_foil_bulk => sumEbb0nuParams,
)
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
    :activity => BkgActivityParams[:Bi214_wire_surface], 
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

singleEParams = Dict(
    :Bi214_foil_bulk => singleEBi214_foil_bulk_Params,
    :Bi214_foil_surface => singleEBi214_foil_surface_Params,
    :Bi214_PMT_glass_bulk => singleEBi214_PMT_glass_bulk_Params,
    :Bi214_wire_bulk => singleEBi214_wire_bulk_Params,
    :Bi214_wire_surface => singleEBi214_wire_surface_Params,
    :Bi210_foil_bulk => singleEBi210_foil_bulk_Params,
    :Bi210_foil_surface => singleEBi210_foil_surface_Params,
    :Bi210_wire_bulk => singleEBi210_wire_bulk_Params,
    :Bi210_wire_surface => singleEBi210_wire_surface_Params,
    :Tl208_foil_bulk => singleETl208_foil_bulk_Params,
    :Tl208_foil_surface => singleETl208_foil_surface_Params,
    :Tl208_PMT_glass_bulk => singleETl208_PMT_glass_bulk_Params,
    :Pa234m_foil_bulk => singleEPa234m_foil_bulk_Params,
    :K40_foil_bulk => singleEK40_foil_bulk_Params,
    :K40_PMT_glass_bulk => singleEK40_PMT_glass_bulk_Params,
    :bb_foil_bulk => singleEbbParams,
    :Xi037_foil_bulk => singleEXiParams,
    :RH037_foil_bulk => singleERH037Params,
    :bb0nu_foil_bulk => singleEbb0nuParams,
)

#### Angular Study


### Bi214
phiBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_bulk], 
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_surface], 
    :bins => binningDict[:Phi],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

phiBi214_PMT_glass_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_glass_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

phiBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

phiBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_surface],
    :bins => binningDict[:Phi],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

### Bi210
phiBi210_foil_bulk_Params = Dict(
    :isotopeName => "Bi210_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_foil_bulk], 
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiBi210_foil_surface_Params = Dict(
    :isotopeName => "Bi210_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_foil_surface], 
    :bins => binningDict[:Phi],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

phiBi210_wire_bulk_Params = Dict(
    :isotopeName => "Bi210_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_wire_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

phiBi210_wire_surface_Params = Dict(
    :isotopeName => "Bi210_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi210_wire_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi210_wire_surface],
    :bins => binningDict[:Phi],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

### Tl208

phiTl208_foil_bulk_Params = Dict(
    :isotopeName => "Tl208_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiTl208_foil_surface_Params = Dict(
    :isotopeName => "Tl208_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_foil_surface],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_surface",
    :amount => SNparams["gasVolume"]
)

phiTl208_PMT_glass_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_glass_bulk], 
    :bins => binningDict[:Phi],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

### Pa234m

phiPa234m_foil_bulk_Params = Dict(
    :isotopeName => "Pa234m_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Pa234m_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### K40

phiK40_foil_bulk_Params = Dict(
    :isotopeName => "K40_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiK40_PMT_glass_bulk_Params = Dict(
    :isotopeName => "K40_PMT_glass_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_glass_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_PMT_glass_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "PMT_glass_bulk",
    :amount => SNparams["PMTGlassMass"]
)

### bb_foil_bulk

phibbParams = Dict(
    :isotopeName => "bb_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### Xi31
phiXiParams = Dict(
    :isotopeName => "Xi037_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Xi037_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### RH037
phiRH037Params = Dict(
    :isotopeName => "RH037_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:RH037_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:RH037_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### bb0nu_foil_bulk

phibb0nuParams = Dict(
    :isotopeName => "bb0nu_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nu_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nu_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiParams = Dict(
    :Bi214_foil_bulk => phiBi214_foil_bulk_Params,
    :Bi214_foil_surface => phiBi214_foil_surface_Params,
    :Bi214_PMT_glass_bulk => phiBi214_PMT_glass_bulk_Params,
    :Bi214_wire_bulk => phiBi214_wire_bulk_Params,
    :Bi214_wire_surface => phiBi214_wire_surface_Params,
    :Bi210_foil_bulk => phiBi210_foil_bulk_Params,
    :Bi210_foil_surface => phiBi210_foil_surface_Params,
    :Bi210_wire_bulk => phiBi210_wire_bulk_Params,
    :Bi210_wire_surface => phiBi210_wire_surface_Params,
    :Tl208_foil_bulk => phiTl208_foil_bulk_Params,
    :Tl208_foil_surface => phiTl208_foil_surface_Params,
    :Tl208_PMT_glass_bulk => phiTl208_PMT_glass_bulk_Params,
    :Pa234m_foil_bulk => phiPa234m_foil_bulk_Params,
    :K40_foil_bulk => phiK40_foil_bulk_Params,
    :K40_PMT_glass_bulk => phiK40_PMT_glass_bulk_Params,
    :bb_foil_bulk => phibbParams,
    :Xi037_foil_bulk => phiXiParams,
    :RH037_foil_bulk => phiRH037Params,
    :bb0nu_foil_bulk => phibb0nuParams,
)