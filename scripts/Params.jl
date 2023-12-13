# include(srcdir("SensitivityModule.jl")) # to get halfLife_to_activity
import SensitivityModule.halfLife_to_activity
#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23,                         # Avogadro's number in [1/mol]
    "W" => 0.08192,                             # Se82 molar mass in [kg/mol]
    "a" => 0.99,                                # abundance/foil enrichment; check number
    "m" => 6.25,                                # foil mass in [kg]
    "t" => 2.5 * 365 * 24 * 3600,               # measurement time in [s]
    "tYear" => 2.5,                             # measurement time in [y]
    "SeThalf" => (9.39 ± 0.17)*1e19 * 365 * 24 * 3600    # Se82 half life in [s], results from NEMO-3
) 

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x
    :Bi214 => (1.50 ± 0.04) / 1000 ,     # [mBq/kg] converted to [Bq/kg]
    :Tl208 => (0.39 ± 0.01) / 1000 ,     # [mBq/kg] converted to [Bq/kg]
    :Pa234m => (17.3 ± 0.1) / 1000 ,    # [mBq/kg] converted to [Bq/kg]
    :K40 => (58.7 ± 0.1) / 1000         # [mBq/kg] converted to [Bq/kg]
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

sumEBi214Params = Dict(
    :isotopeName => "Bi214", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumETl208Params = Dict(
    :isotopeName => "Tl208", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEPa234mParams = Dict(
    :isotopeName => "Pa234m", 
    :signal => :false, 
    :activity => BkgActivityParams[:Pa234m], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEK40Params = Dict(
    :isotopeName => "K40", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

sumEbbParams = Dict(
    :isotopeName => "2nubb", 
    :signal => :true, 
    :activity => SigActivityParams[:bb2Standard], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => sumEParams[:nTotalSim], 
    :bins => sumEParams[:binning]
)

#### Single energy spectra
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