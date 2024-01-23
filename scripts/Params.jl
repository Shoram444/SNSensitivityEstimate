# include(srcdir("SensitivityModule.jl")) # to get halfLife_to_activity
import SensitivityModule.halfLife_to_activity
using Measurements
#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23,                         # Avogadro's number in [1/mol]
    "W" => 0.08192,                             # Se82 molar mass in [kg/mol]
    "a" => 0.99,                                # abundance/foil enrichment; check number
    "foilMass" => 6.25,                         # foil mass in [kg]
    "gasVolume" => 6*2*4,                       # tracker volume in [m3]
    "PMTGlassMass" => 286,                      # PMT glass mass in [kg]
    "wireBulkMass" => 20,                       # DUMMY VALUE!! mass of the tracker wires in [kg]
    "t" => 2.5 * 365 * 24 * 3600,               # measurement time in [s]
    "tYear" => 2.5,                             # measurement time in [y]
    "SeThalf" => (9.39 ± 0.17)*1e19 * 365 * 24 * 3600    # Se82 half life in [s], results from NEMO-3
) 

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x and docDB 4505
    :Bi214_foil_bulk => 300 / 1_000_000 ,       # [μBq/kg] converted to [Bq/kg] from SN measurements
    :Tl208_foil_bulk => 54 / 1_000_000 ,        # [μBq/kg] converted to [Bq/kg] from SN measurements
    :Pa234m_foil_bulk => (17.3 ± 0.1) / 1000 ,  # [mBq/kg] converted to [Bq/kg] from NEMO-3
    :K40_foil_bulk => (58.7 ± 0.1) / 1000,      # [mBq/kg] converted to [Bq/kg] from NEMO-3
    :Bi214_PMT_bulk => 140 / 286 ,              # [Bq/kg] originally the value is given as 140Bq, I just divide by PMT weight here
    :Tl208_PMT_bulk => 41.4 / 286 ,             # [Bq/kg] originally the value is given as 41.4Bq, I just divide by PMT weight here
    :K40_PMT_bulk => 417 / 286,                 # [Bq/kg] originally the value is given as 417Bq, I just divide by PMT weight here
    :Bi214_wire_surface => 150 / 1_000_000 ,           # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3) 
    :Bi214_wire_bulk => 0.00001                 # MOCK VALUE
)

SigActivityParams = Dict( 
    :bb2Standard => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf"]), # activity calculated from 2nubb Half-life in [Bq/kg]
    :Xi037 => (1.5 ± 0.1) / 1000 # mock value
)

SimulationParams = Dict(
    :Bi214_foil_bulk => 1e8,
    :Bi214_foil_surface => 1e7,
    :Bi214_PMT_bulk => 1e7,
    :Bi214_wire_bulk => 1e7,
    :Bi214_wire_surface => 1e7,

    :Tl208_foil_bulk => 1e8,
    :Tl208_PMT_bulk => 1e8,

    :Pa234m_foil_bulk => 1e8,

    :K40_foil_bulk => 1e8,
    :K40_PMT_bulk => 1e8,

    :bb2Stabdard_foil_bulk => 1e8,

    :Xi037_foil_bulk => 1e8 
)

#### Sum energy spectra
sumEParams = Dict(
    :binning => 0:100:3500, 
)

### Bi214
sumEBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_bulk], 
    :bins => sumEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumEBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_surface], 
    :bins => sumEParams[:binning],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]

)

sumEBi214_PMT_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_bulk],
    :bins => sumEParams[:binning],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

sumEBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_bulk],
    :bins => sumEParams[:binning],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

sumEBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_surface],
    :bins => sumEParams[:binning],
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
    :bins => sumEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumETl208_PMT_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_bulk], 
    :bins => sumEParams[:binning],
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
    :bins => sumEParams[:binning],
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
    :bins => sumEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

sumEK40_PMT_bulk_Params = Dict(
    :isotopeName => "K40_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_PMT_bulk],
    :bins => sumEParams[:binning],
    :vertexPosition => "PMT_bulk",
    :amount => SNparams["PMTGlassMass"]
)

### 2nubb

sumEbbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb2Stabdard_foil_bulk],
    :bins => sumEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### Xi31
sumEXiParams = Dict(
    :isotopeName => "Xi31", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Xi037_foil_bulk],
    :bins => sumEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)


###############################################################################
###############################################################################
#### Single energy spectra
###############################################################################
###############################################################################

singleEParams = Dict(
    :binning => 0:100:3500,
    :nTotalSim => 2e8           # because we have 2 electrons per event
)

### Bi214
singleEBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_bulk], 
    :bins => singleEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleEBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_surface], 
    :bins => singleEParams[:binning],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

singleEBi214_PMT_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_bulk],
    :bins => singleEParams[:binning],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

singleEBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_bulk],
    :bins => singleEParams[:binning],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

singleEBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_surface],
    :bins => singleEParams[:binning],
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
    :bins => singleEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleETl208_PMT_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_bulk], 
    :bins => singleEParams[:binning],
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
    :bins => singleEParams[:binning],
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
    :bins => singleEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

singleEK40_PMT_bulk_Params = Dict(
    :isotopeName => "K40_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_PMT_bulk],
    :bins => singleEParams[:binning],
    :vertexPosition => "PMT_bulk",
    :amount => SNparams["PMTGlassMass"]
)

### 2nubb

singleEbbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb2Stabdard_foil_bulk],
    :bins => singleEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### Xi31
singleEXiParams = Dict(
    :isotopeName => "Xi31", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Xi037_foil_bulk],
    :bins => singleEParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)


#### Angular Study
phiParams = Dict(
    :binning => 0:5:180,
    :nTotalSim => 1e8
)

### Bi214
phiBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_bulk], 
    :bins => phiParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_foil_surface], 
    :bins => phiParams[:binning],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

phiBi214_PMT_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_bulk],
    :bins => phiParams[:binning],
    :vertexPosition => "PMT_glass",
    :amount => SNparams["PMTGlassMass"]
)

phiBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_bulk],
    :bins => phiParams[:binning],
    :vertexPosition => "wire_bulk",
    :amount => SNparams["wireBulkMass"]
)

phiBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_wire_surface],
    :bins => phiParams[:binning],
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
    :bins => phiParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiTl208_PMT_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_PMT_bulk], 
    :bins => phiParams[:binning],
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
    :bins => phiParams[:binning],
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
    :bins => phiParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

phiK40_PMT_bulk_Params = Dict(
    :isotopeName => "K40_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_PMT_bulk],
    :bins => phiParams[:binning],
    :vertexPosition => "PMT_bulk",
    :amount => SNparams["PMTGlassMass"]
)

### 2nubb

phibbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb2Stabdard_foil_bulk],
    :bins => phiParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### Xi31
phiXiParams = Dict(
    :isotopeName => "Xi31", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Xi037_foil_bulk],
    :bins => phiParams[:binning],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)
