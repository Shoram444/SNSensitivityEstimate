using DrWatson
@quickactivate "SNSensitivityEstimate"

println("loading pkgs")

lprob_to_prob(x) = abs(10^(-x))
prob_to_lprob(x) = abs(log10(x))

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Metaheuristics

println("loaded pkgs")


## Helper functions
include(scriptsdir("phd_final/07_sensitivity/sensitivity_nd/helper_functions.jl"))


analysisDict = Dict(
    :signal => "Nnubb500keV_foil_bulk",
    :radon_tag => 1,
    :hours => 22
)

const VARS_ = [
    "phi", 
    "sumE", 
    # "maxE", 
    # "minE",
    # "avgE", 
    # "r", 
    # "singleE",
    "dy", 
    "dz",
    # "sameSide",
    # "Pint",
    # "Pext",
    "lPint",
    "lPext",
    # "trackLength1",
    # "trackLength2",
    # "caloTime1",
    # "caloTime2",
    ]

const BOUNDS_ = (
    phi = (0,180),
    sumE = (300, 3000),
    # minE = (0, 3500),
    # maxE = (0, 3500),
    # avgE = (0, 3500),
    # r = (0, 200),
    # singleE = (0, 3500),
    dy = (0, 150),
    dz = (0, 150),
    # sameSide = (0, 1),
    # Pint = (0, 1),
    # Pext = (0, 1),
    lPint = (0, 10),
    lPext = (0, 50),
    # trackLength1 = (0, 3000),
    # trackLength2 = (0, 3000),
    # caloTime1 = (0, 100),
    # caloTime2 = (0, 100),
)

processes = load_ndim_processes(datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41"), BOUNDS_, VARS_)

data_dir = datadir("sims/final_phd/fal5_12perc_Boff_Cimrman_J41/neutrons_jan_2026/")

include(joinpath(data_dir, "read_neutrons_ND.jl"))
neutron_processes = load_neutron_processes_ND(data_dir, VARS_, BOUNDS_)
for p in neutron_processes
    set_activity!(p, p.activity / 3)
    println("Neutron process: ", p.isotopeName, " nTotalSim: ", p.nTotalSim, " activity: ", p.activity)
end


# signal_name = "Nnubb500keV_foil_bulk"
signal_name = analysisDict[:signal]
signal = get_process(signal_name, processes) |> first
# signal = get_process("RH037_foil_bulk", processes) |> first
# signal = get_process("bb0nuM2_foil_bulk", processes)

# declare background processes
bkgs = [    
    # 2nubb
    "bb_foil_bulk", 

    # internal contamination
    "Bi214_foil_bulk", 
    "K40_foil_bulk", 
    "Pa234m_foil_bulk", 
    "Tl208_foil_bulk", 

    # radon
    "Bi210_wire_surface",     
    "Bi214_wire_surface", 

    # detector contamination
    "Bi214_PMT_glass_bulk",
    "Bi214_PMT5_glass_bulk",
    "K40_PMT_glass_bulk", 
    "K40_PMT5_glass_bulk", 
    "Tl208_PMT_glass_bulk", 
    "Tl208_PMT5_glass_bulk", 
    "Bi214_calo_8inch_scin_bulk", 
    "K40_calo_8inch_scin_bulk", 

    # gamma
    "Tl208_hall_bulk", 
    "Bi214_hall_bulk", 
    "K40_hall_bulk",
]


background = [get_process(bkg, processes) |> first for bkg in bkgs]
background = vcat(background, neutron_processes)

# set 2nubb to background process (initially it's signal for exotic 2nubb analyses)
set_signal!(background[1], false)
# set_activity!(background[7], 0.15/1e3) # radon to 0.15mBq/kg
# set_activity!(background[7], 2/1e3) # radon to 2mBq/kg
radon_tag = analysisDict[:radon_tag]
if radon_tag == 1
    set_activity!(background[7], 150/1e6) # radon to 150 uBq/kg
elseif radon_tag == 2
    set_activity!(background[7], 2/1e3) # radon to 2 mBq/kg
elseif radon_tag == 3
    set_activity!(background[7], 0.6/1e3) # radon to 0.6 mBq/kg
end

const α_ = 1.64485362695147

println("loaded files, signal = $(signal.isotopeName)")


function roi_vector_to_normalized(roi::AbstractVector)

    x̂ = Vector{Float64}(undef, length(roi))
    WMIN = 0.05   # minimum width fraction

    for (i, var) in enumerate(VARS_)
        L, U = BOUNDS_[Symbol(var)]
        R = U - L

        minv = roi[2i-1]
        maxv = roi[2i]

        s = (minv - L) / R
        w = (maxv - minv) / (U - minv)

        x̂[2i-1] = clamp(s, 0.0, 1.0)
        x̂[2i]   = clamp((w - WMIN) / (1 - WMIN), 0.0, 1.0)
    end

    return x̂
end

function normalized_to_roi_vector(x̂::AbstractVector)
    roi = Vector{Float64}(undef, length(x̂))
    WMIN = 0.05   # minimum width fraction

    for (i, var) in enumerate(VARS_)
        L, U = BOUNDS_[Symbol(var)]
        R = U - L

        ŝ = clamp(x̂[2i-1], 0.0, 1.0)   # start fraction
        ŵ = clamp(x̂[2i],   0.0, 1.0)   # width fraction

        w = WMIN + (1 - WMIN) * ŵ       # enforce minimum width

        minv = L + ŝ * R
        maxv = minv + w * (U - minv)

        roi[2i-1] = minv
        roi[2i]   = maxv
    end

    return roi
end



searchRange = make_stepRange(signal)

# n = 100
# lb = first.(searchRange)
# ub = last.(searchRange)

# manual lower upper bounds
lb = float.([0, 140, 300, 2000, 0, 50, 0, 50, 0, 0, 0, 50 ])
ub = float.([80, 180, 1600, 2800, 10, 150, 10, 150, 0, 10, 10, 50])


const PROCESSES = vcat(signal, background)
const SNparams_ = SNparams

function prob(x::Vector{Float64})
    # xs = normalized_to_roi_vector(collect(x))
    return - get_s_to_b(SNparams_, α_, PROCESSES, x; approximate="formula")
end



options = Options(;
    x_tol = 1e-2,
    f_tol = 1e-4,
    f_tol_rel = 1e-5,
    f_tol_abs = 1e-5,
    time_limit = analysisDict[:hours]*60*60.0, #60*60*12.0,
    parallel_evaluation = false,
    verbose = true,
)


bounds = boxconstraints(lb = lb, ub = ub)

algo = ECA(N=20, K=3, ;options)

x0 = ( float.( [0, 180, 400, 2400, 0, 120, 0, 120, 0, 4, 1, 50] ) )
set_user_solutions!(algo, x0, prob)


result = Metaheuristics.optimize(prob, bounds, algo)
@show minimum(result)
@show res=  minimizer(result)

best = get_best_ROI_ND(res, signal)
best_sens = get_sensitivityND(SNparams, α, vcat(signal, background), best; approximate="table")

print(best_sens)


