import SNSensitivityEstimate.halfLife_to_activity
using Measurements

include("SNparams.jl") 
include("Activities.jl") 
include("NSimulated.jl") 



binningDict = Dict(
    :SumE => 0:100:3500,
    :SingleE => 0:100:3500,
    :Phi => 0:5:180,
)

include("SumE.jl") 
include("SingleE.jl") 
include("Phi.jl") 


