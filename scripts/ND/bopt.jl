using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

push!(LOAD_PATH, srcdir())
using SensitivityModule, CairoMakie, Surrogates 


# File "scripts/Params.jl" contains the all (most) of the necessary parameters for the sensitivity estimation in one place
# Information is placed in `Dict` (Dictionaries). Take a look inside for details, but the general idea is we export these 
# dictionaries into this script, which uses their values. 
println("loaded pkgs")
include(scriptsdir("Params.jl"))

vars = [
    "phi", 
    "sumE", 
    # "maxE", 
    # "minE", 
    "r", 
    # "dy", 
    # "dz",
    # "sameSide"
    ]

bins = (
    phi = (0,180),
    sumE = (0, 3500),
    # maxE = (0, 3500),
    # minE = (0, 3500),
    r = (0, 100),
    # dy = (-100, 100),
    # dz = (-100, 100),
    # sameSide = (0, 1)
)

processes = load_ndim_processes("fal5_TKrec", bins, vars)

signal = get_process("bb0nu_foil_bulk", processes)
# signal = get_process("RH037_foil_bulk", processes)
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
background = [
    get_process("bb_foil_bulk", processes),
    get_process("Bi214_foil_bulk", processes),
    get_process("Bi214_wire_surface", processes),
    get_process("Tl208_foil_bulk", processes),
    get_process("K40_foil_bulk", processes),
    get_process("Pa234m_foil_bulk", processes),
]

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)

# set the number of total simulated events (there's a default in "scripts/Params.jl", but this is usecase dependend)
set_nTotalSim!( signal, 1e8 )
# set_nTotalSim!( signal, 1e8 )
set_nTotalSim!( background[1], 0.99e8 )
set_nTotalSim!( background[2], 1e8 )
set_nTotalSim!( background[3], 1e8 )
set_nTotalSim!( background[4], 1e8 )
set_nTotalSim!( background[5], 1e8 )
set_nTotalSim!( background[6], 1e8 )

α= 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


prob(x) = - SensitivityModule.get_s_to_b(SNparams, α, vcat(signal, background), x; approximate="formula")


function make_stepRange(process)
    stepRange = Tuple{Int64, Int64}[]
    for k in keys(process.bins) 
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
        push!(stepRange, (process.bins[k][1], process.bins[k][2]))
    end
    return stepRange
end


searchRange = make_stepRange(signal)

n_samples = 6000
lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float

# lower_bound = [
#     0, 160,
#     0, 1000,
#     0, 1000,
#     0, 30
# ] .|> float

# upper_bound = [
#     30, 180,
#     1000, 3500,
#     1000, 3500,
#     30, 50
# ] .|> float


xys = Surrogates.sample(n_samples, lower_bound, upper_bound, SobolSample())
xs = [[x[i] for i in 1:length(lower_bound)] for x in xys]
zs = [prob(x) for x in xs]
replace!(zs, 1e20 => 0.0)
zs = Measurements.value.(zs)

phi1 = [x[1] for x in xs]
phi2 = [x[2] for x in xs]
scatter(phi1, phi2, zs, axis=(type=Axis3, xlabel = "phi1", ylabel = "phi2", zlabel = "s/b"))

esum1 = [x[3] for x in xs]
esum2 = [x[4] for x in xs]
scatter(esum1, esum2, zs, axis=(type=Axis3, xlabel = "esum1", ylabel = "esum2", zlabel = "s/b"))

# emin1 = [x[5] for x in xs]
# emin2 = [x[6] for x in xs]
# scatter(emin1, emin2, zs, axis=(type=Axis3, xlabel = "emin1", ylabel = "emin2", zlabel = "s/b"))


r1 = [x[5] for x in xs]
r2 = [x[6] for x in xs]
scatter(r1, r2, zs, axis=(type=Axis3, xlabel = "r1", ylabel = "r2", zlabel = "s/b"))


using AbstractGPs, SurrogatesAbstractGPs
krig = Kriging(xys, zs, lower_bound, upper_bound, p = fill(2.0, 8), theta = fill(0.03, 8))
result = Surrogates.surrogate_optimize(prob, EI(), lower_bound, upper_bound, krig, SobolSample(); maxiters = 100)
res = collect(result[1])

function get_best_ROI_ND(res::Vector{<:Real}, process)
    best = res
    best_roi = NamedTuple(
        k => (round(best[i]), round(best[i+1])) 
        for (i,k) in zip(1:2:length(process.bins)*2-1, keys(process.bins))
    )
    return best_roi
end

# best = get_best_ROI_ND(res, signal)
# get_sensitivityND(SNparams, α, vcat(signal, background), best; approximate="formula")

let 
    fig = Figure(size = (800, 1200))
    ### phi
    ax1 = Axis(fig[1,1], title="phi1 vs phi2", xlabel="phi1", ylabel="phi2")
    scatter!(ax1, phi1, phi2, color=zs, alpha= 0.7, markersize = 6)

    ax3 = Axis(fig[1,2], title="phi1 vs z", xlabel="phi1", ylabel="z")
    scatter!(ax3, phi1, zs)

    ax3 = Axis(fig[2,1], title="phi2 vs z", xlabel="phi2", ylabel="z")
    scatter!(ax3, phi2, zs)

    ### esum
    ax2 = Axis(fig[3,1], title="esum1 vs esum2", xlabel="esum1", ylabel="esum2")
    scatter!(ax2, esum1, esum2, color=zs, alpha= 0.7, markersize = 6)

    ax4 = Axis(fig[3,2], title="esum1 vs z", xlabel="esum1", ylabel="z")
    scatter!(ax4, esum1, zs)

    ax4 = Axis(fig[4,1], title="esum2 vs z", xlabel="esum2", ylabel="z")
    scatter!(ax4, esum2, zs)

    ### r
    ax2 = Axis(fig[5,1], title="r1 vs r2", xlabel="r1", ylabel="r2")
    scatter!(ax2, r1, r2, color=zs, alpha= 0.7, markersize = 6)

    ax4 = Axis(fig[5,2], title="r1 vs z", xlabel="r1", ylabel="z")
    scatter!(ax4, r1, zs)

    ax4 = Axis(fig[6,1], title="r2 vs z", xlabel="r2", ylabel="z")
    scatter!(ax4, r2, zs)
    save( "plots/heuristic_3d_plots.png",fig, px_per_unit = 6)
    fig
    
end
upper_bound
let 
    e1, e2 = range(0, 3500, 60), range(0, 3500, 60)
    m = [e1 e2]
    k(m) = krig([0.0 180.0 m[1] m[2] 0.0 100.0] )

    fig = Figure(size = (800, 600))
    ax1 = Axis3(fig[1,1], azimuth = 0.3*pi)

    surface!(ax1, e1, e2, (e1,e2) -> k([e1 e2]))
    fig
end

let 
    phi1, phi2 = range(0, 180, 20), range(0, 180, 20)
    m = [phi1 phi2]
    k(m) = krig([m[1] m[2] 0.0 3500.0 0.0 100.0])

    fig = Figure(size = (800, 600))
    ax1 = Axis3(fig[1,1], azimuth = 0.1*pi)

    surface!(ax1, phi1, phi2, (phi1,phi2) -> k([phi1 phi2]))
    fig
end

e1, e2 = range(0, 3500, 20), range(0, 3500, 20)
m = [e1 e2]
k(m) = krig([0.0 180.0 m[1] m[2]])

k([e1[2] e2[2]])









