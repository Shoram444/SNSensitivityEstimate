# include(srcdir("SensitivityModule.jl")) # to get halfLife_to_activity
import SensitivityModule.halfLife_to_activity
using Measurements
#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23,                         # Avogadro's number in [1/mol]
    "W" => 0.08192,                             # Se82 molar mass in [kg/mol]
    "a" => 0.99,                                # abundance/foil enrichment; check number
    "m" => 6.25,                                # foil mass in [kg]
    "V" => 6*2*4,                               # tracker volume in [m3]
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
    :Bi214_radon => 150 / 1_000_000 ,           # [μBq/m3] converted to [Bq/m3] from SN measurements; upper limit 48m3 (4x6x2m3) 
)

SigActivityParams = Dict( 
    :bb2Standard => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf"]), # activity calculated from 2nubb Half-life in [Bq/kg]
    :Xi037 => (1.5 ± 0.1) / 1000 # mock value
)

#### Sum energy spectra
sumEParams = Dict(
    :binning => 0:100:3500, 
    :nTotalSim => 1e8
)

### Bi214
sumEBi214_foil_bulk_Params = Dict(
    :isotopeName => "Bi214_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEBi214_foil_surface_Params = Dict(
    :isotopeName => "Bi214_foil_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEBi214_PMT_bulk_Params = Dict(
    :isotopeName => "Bi214_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEBi214_wire_bulk_Params = Dict(
    :isotopeName => "Bi214_wire_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEBi214_wire_surface_Params = Dict(
    :isotopeName => "Bi214_wire_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_radon], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

### Tl208

sumETl208_foil_bulk_Params = Dict(
    :isotopeName => "Tl208_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumETl208_PMT_bulk_Params = Dict(
    :isotopeName => "Tl208_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

### Pa234m

sumEPa234m_foil_bulk_Params = Dict(
    :isotopeName => "Pa234m_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

### K40

sumEK40_foil_bulk_Params = Dict(
    :isotopeName => "K40_foil_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEK40_PMT_bulk_Params = Dict(
    :isotopeName => "K40_PMT_bulk", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_PMT_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

### 2nubb

sumEbbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

### Xi31
sumEXiParams = Dict(
    :isotopeName => "Xi31", 
    :signal => :true, 
    :activity => SigActivityParams[:Xi037], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
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

singleEBi214Params = Dict(
    :isotopeName => "Bi214", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => singleEParams[:nTotalSim], 
    :bins => singleEParams[:binning]
)

singleETl208Params = Dict(
    :isotopeName => "Tl208", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => singleEParams[:nTotalSim], 
    :bins => singleEParams[:binning]
)

singleEPa234mParams = Dict(
    :isotopeName => "Pa234m", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => singleEParams[:nTotalSim], 
    :bins => singleEParams[:binning]
)

singleEK40Params = Dict(
    :isotopeName => "K40", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => singleEParams[:nTotalSim], 
    :bins => singleEParams[:binning]
)

singleEbbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => singleEParams[:nTotalSim], 
    :bins => singleEParams[:binning]
)

#### Angular Study
phiParams = Dict(
    :binning => 0:5:180,
    :nTotalSim => 1e8
)

phiBi214Params = Dict(
    :isotopeName => "Bi214", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => phiParams[:nTotalSim], 
    :bins => phiParams[:binning]
)

phiTl208Params = Dict(
    :isotopeName => "Tl208", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => phiParams[:nTotalSim], 
    :bins => phiParams[:binning]
)

phiPa234mParams = Dict(
    :isotopeName => "Pa234m", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => phiParams[:nTotalSim], 
    :bins => phiParams[:binning]
)

phiK40Params = Dict(
    :isotopeName => "K40", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => phiParams[:nTotalSim], 
    :bins => phiParams[:binning]
)

phibbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => phiParams[:nTotalSim], 
    :bins => phiParams[:binning]
)