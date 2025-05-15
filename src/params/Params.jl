import SensitivityModule.halfLife_to_activity
using Measurements

include(srcdir("params/SNparams.jl")) 
include(srcdir("params/Activities.jl")) 
include(srcdir("params/NSimulated.jl")) 



binningDict = Dict(
    :SumE => 0:100:3500,
    :SingleE => 0:100:3500,
    :Phi => 0:5:180,
)

include(srcdir("params/SumE.jl")) 
include(srcdir("params/SingleE.jl")) 
include(srcdir("params/Phi.jl")) 


