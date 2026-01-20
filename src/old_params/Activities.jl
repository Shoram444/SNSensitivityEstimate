# Background parameters

BkgActivityParams = Dict( #activities from Table 1 from 10.1140/epjc/s10052-018-6295-x and docDB 4505
    :Pa234m_foil_bulk => 17.3 / 1000 ,          # [Bq/kg] NEMO3
    :Bi214_foil_bulk => 10 / 1_000_000,               # [Bq/kg] from SN target of Bi214, assume same
    :Bi214_foil_surface => 150 / 1_000_000 ,                  # [Bq/m3] from SN target; gas volume 15m3 
    :Bi214_wire_surface => 150 / 1_000_000 ,                  # [Bq/m3] from SN target; gas volume 15m3 
    :Bi214_field_wires => 150 / 1_000_000 ,                  # [Bq/m3] from SN target; gas volume 15m3 
    :Bi214_wire_bulk => 0.00001 ,                       # MOCK VALUE
    :Bi214_PMT_glass_bulk => 141 / 286,        # [Bq/kg] originally the value is given as 141Bq, I just divide by PMT weight here, from Malak thesis
    :Tl208_foil_bulk => 2 / 1_000_000,                # [Bq/kg] from SN target
    :Tl208_foil_surface => 150 / 1_000_000 ,                  # [Bq/m3] from SN target; gas volume 15m3
    :Tl208_PMT_glass_bulk => 115 / 286 ,               # [Bq/kg] originally the value is given as 115Bq, I just divide by PMT weight here, from Malak thesis
    :Bi210_foil_bulk => 10 / 1_000_000 ,                # MOCK VALUE
    :Bi210_foil_surface => 150 / 1_000_000 ,                  # [Bq/m3] from SN target; gas volume 15m3 
    :Bi210_wire_surface => 150 / 1_000_000 ,                  # [Bq/m3] from SN target; gas volume 15m3 
    :Bi210_wire_bulk => 0.00001 ,                       # MOCK VALUE
    :K40_foil_bulk => 58.7 / 1000 ,             # [Bq/kg] NEMO3
    :K40_PMT_glass_bulk => 417 / 286 ,                  # [Bq/kg] originally the value is given as 417Bq, I just divide by PMT weight here, from Malak thesis
    :Bi214_hall_surface => 0.5 ,                  # [Bq] Value is taken from Xalbat (flux) divided twice by 500 - for tof cut and for shielding effect 
    :Tl208_hall_surface => 2.4 ,                  # [Bq] Value is taken from Xalbat (flux) divided twice by 500 - for tof cut and for shielding effect 
    :K40_hall_surface => 1.2 ,                  # [Bq] Value is taken from Xalbat (flux) divided twice by 500 - for tof cut and for shielding effect 
    :gamma_experimental_surface => 1 / 430,         # mock gammas for "somewhat flat" spectrum in nu0 ROI
    :gamma_hall_bulk => 1 / 430,         # mock gammas for "somewhat flat" spectrum in nu0 ROI
    :K40_calo_8inch_scin_bulk => 2.2/1000,         # [Bq/kg] from calo Federic measurements internal, docdb 3174-v2
    :Bi214_calo_8inch_scin_bulk => 0.3/1000,         # [Bq/kg] from calo Federic measurements internal, docdb 3174-v2
    )


# Signal parameters
SigActivityParams = Dict( 
    :bb_foil_bulk =>        halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # activity calculated from 2nubb Half-life in [Bq/kg]
    :Xi037_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH037_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH020_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH025_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH030_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH035_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH040_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH045_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :RH050_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :bb0nuRHl_foil_bulk =>  halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"]), # mock value
    :bb0nu_foil_bulk =>     halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
    :bb0nuM1_foil_bulk =>   halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
    :bb0nuM2_foil_bulk =>   halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
    :bb0nuScalar0_foil_bulk =>   halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
    :bb0nuScalar1_foil_bulk =>   halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
    :bb0nuScalar2_foil_bulk =>   halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf0nu"]), # mock value
)







