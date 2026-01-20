#### SN sensitivity parameters ####
SNparams = Dict(
    "Nₐ" => 6.02214e23,                         # Avogadro's number in [1/mol]
    "W" => 0.08192,                             # Se82 molar mass in [kg/mol]
    "a" => 0.9776,                               # abundance/foil enrichment; check number
    "foilMass" => 6.25,                         # foil mass in [kg]
    "gasVolume" => 15,                          # tracker volume in [m3]
    "PMTGlassMass" => 286,                      # PMT glass mass in [kg]
    "wireBulkMass" => 20,                       # DUMMY VALUE!! mass of the tracker wires in [kg]
    "hall_surface" => 1,                       # External gamma activity is given in [Bq] so no need to scale by amount
    "hall_bulk" => 1,                       # External gamma activity is given in [Bq] so no need to scale by amount
    "t" => 2.86 * 365 * 24 * 3600,              # measurement time in [s]
    "tYear" => 2.86,                             # measurement time in [y]
    "SeThalf2nu" => 9.39 * 1e19 * 365 * 24 * 3600,    # 2nu Se82 half life in [s], results from NEMO-3
    "SeThalf0nu" => 1e26 * 365 * 24 * 3600,     # 0nu Se82 half life in [s], results from NEMO-3
    "Q" => 2997.9,                               # [keV]
    "caloMass_8inch" => 10.0,                    # mass of one 8 inch calo scintillator in [kg]
) 