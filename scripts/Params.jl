include(srcdir("SensitivityModule.jl")) # to get halfLife_to_activity
#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23,                         # Avogadro's number in [1/mol]
    "W" => 0.08192,                             # Se82 molar mass in [kg/mol]
    "a" => 0.99,                                # abundance/foil enrichment; check number
    "m" => 6.25,                                # foil mass in [kg]
    "t" => 2.5 * 365 * 24 * 3600,               # measurement time in [s]
    "SeThalf" => 9.39*1e19 * 365 * 24 * 3600    # Se82 half life in [s], results from NEMO-3
) 

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x
    :Bi214 => 1.50 / 1000 ,     # [mBq/kg] converted to [Bq/kg]
    :Tl208 => 0.39 / 1000 ,     # [mBq/kg] converted to [Bq/kg]
    :Pa234m => 17.3 / 1000 ,    # [mBq/kg] converted to [Bq/kg]
    :K40 => 58.7 / 1000         # [mBq/kg] converted to [Bq/kg]
)

SigActivityParams = Dict( 
    :bb2Standard => halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf"]), # activity calculated from 2nubb Half-life in [Bq/kg]
)

#### Sum energy spectra
sumEParams = Dict(
    :binning => 0:100:3500, 
    :nTotalSim => 1e8
)

#### Single energy spectra
singleEParams = Dict(
    :binning => 0:100:3500,
    :nTotalSim => 2e8           # because we have 2 electrons per event
)

#### Angular Study
phiParams = Dict(
    :binning => 0:3:180,
    :nTotalSim => 1e8
)
