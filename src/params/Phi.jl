
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

phiBi214_field_wires_Params = Dict(
    :isotopeName => "Bi214_field_wires", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_field_wires], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_field_wires],
    :bins => binningDict[:Phi],
    :vertexPosition => "radon",
    :amount => SNparams["gasVolume"]
)

phiBi214_hall_surface_Params = Dict(
    :isotopeName => "Bi214_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Bi214_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Bi214_hall_surface],
    :bins => binningDict[:Phi],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
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

phiTl208_hall_surface_Params = Dict(
    :isotopeName => "Tl208_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:Tl208_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:Tl208_hall_surface],
    :bins => binningDict[:Phi],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
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

phiK40_hall_surface_Params = Dict(
    :isotopeName => "K40_hall_surface", 
    :signal => :false, 
    :activity => BkgActivityParams[:K40_hall_surface], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:K40_hall_surface],
    :bins => binningDict[:Phi],
    :vertexPosition => "hall_surface",
    :amount => SNparams["hall_surface"]
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

### bb0nuM1_foil_bulk

phibb0nuM1Params = Dict(
    :isotopeName => "bb0nuM1_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuM1_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuM1_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### bb0nu_foil_bulk

phibb0nuM2Params = Dict(
    :isotopeName => "bb0nuM2_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuM2_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuM2_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### bb0nuRHl_foil_bulk

phibb0nuRHlParams = Dict(
    :isotopeName => "bb0nuRHl_foil_bulk", 
    :signal => :true, 
    :activity => SigActivityParams[:bb0nuRHl_foil_bulk], 
    :timeMeas => SNparams["t"], 
    :nTotalSim => SimulationParams[:bb0nuRHl_foil_bulk],
    :bins => binningDict[:Phi],
    :vertexPosition => "foil_bulk",
    :amount => SNparams["foilMass"]
)

### neutron_external
phiNeutronParams = Dict(
    :isotopeName => "neutron_external",
    :signal => :false,
    :bins => binningDict[:Phi],
)

phiParams = Dict(
    :Bi214_foil_bulk => phiBi214_foil_bulk_Params,
    :Bi214_foil_surface => phiBi214_foil_surface_Params,
    :Bi214_hall_surface => phiBi214_hall_surface_Params,
    :Bi214_PMT_glass_bulk => phiBi214_PMT_glass_bulk_Params,
    :Bi214_wire_bulk => phiBi214_wire_bulk_Params,
    :Bi214_wire_surface => phiBi214_wire_surface_Params,
    :Bi214_field_wires => phiBi214_field_wires_Params,
    :Bi210_foil_bulk => phiBi210_foil_bulk_Params,
    :Bi210_foil_surface => phiBi210_foil_surface_Params,
    :Bi210_wire_bulk => phiBi210_wire_bulk_Params,
    :Bi210_wire_surface => phiBi210_wire_surface_Params,
    :Tl208_foil_bulk => phiTl208_foil_bulk_Params,
    :Tl208_foil_surface => phiTl208_foil_surface_Params,
    :Tl208_hall_surface => phiTl208_hall_surface_Params,
    :Tl208_PMT_glass_bulk => phiTl208_PMT_glass_bulk_Params,
    :Pa234m_foil_bulk => phiPa234m_foil_bulk_Params,
    :K40_foil_bulk => phiK40_foil_bulk_Params,
    :K40_hall_surface => phiK40_hall_surface_Params,
    :K40_PMT_glass_bulk => phiK40_PMT_glass_bulk_Params,
    :bb_foil_bulk => phibbParams,
    :Xi037_foil_bulk => phiXiParams,
    :RH037_foil_bulk => phiRH037Params,
    :bb0nu_foil_bulk => phibb0nuParams,
    :bb0nuM1_foil_bulk => phibb0nuM1Params,
    :bb0nuM2_foil_bulk => phibb0nuM2Params,
    :bb0nuRHl_foil_bulk => phibb0nuRHlParams,
    :neutron_external => phiNeutronParams
)