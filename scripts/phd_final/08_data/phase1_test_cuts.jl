using DrWatson
@quickactivate "SNSensitivityEstimate"

using Revise
using SNSensitivityEstimate, CairoMakie, DataFramesMeta, CSV, Random, FHist, Distributions
using UnROOT

println("loaded pkgs")

include(datadir("sims/final_phd/data_sims/load_nd_sim_files.jl"))


function filter_data!(d, roi)
    d_new = deepcopy(d)
    for key in keys(roi)
        lo, hi = roi[key]
        filter!(row -> lo <= row[key] <= hi, d_new)
    end
    return d_new
end

function smear_energy(E::Real, fwhm::Real) 
    if(fwhm == 0)
        return E
    end
    sigma = fwhm/2.355 * E * sqrt(1000/E)
    return rand(Normal(E, sigma))
end

function get_smeared_energy(p, var, fwhm)
    if var == :simulatedEnergy1 || var == :simulatedEnergy2
        return smear_energy.(getproperty.(p.data, var), fwhm)
    elseif var == :sumEsimu
        return smear_energy.(getproperty.(p.data, :simulatedEnergy1), fwhm) + smear_energy.(getproperty.(p.data, :simulatedEnergy2), fwhm)
    else
        return getproperty.(p.data, var)
    end
end

function get_expected_simu_counts(p, var, fwhm; bins = (0:100:4000))
    smeared_energy = fwhm == 0 ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm)
    h = normalize(Hist1D(smeared_energy; binedges = bins), width = false)
    n_exp = p.activity * p.timeMeas * p.amount * (length(p.data) / p.nTotalSim)
    return h .* n_exp
end


f1 = ROOTFile(datadir("data/final_phd/phase1_data.root"))
d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

f2 = ROOTFile(datadir("data/final_phd/phase2_data.root"))
d2 = LazyTree(f2, "tree", keys(f2["tree"])) |> DataFrame

f3 = ROOTFile(datadir("data/final_phd/phase3_data.root"))
d3 = LazyTree(f3, "tree", keys(f3["tree"])) |> DataFrame

begin

    data_cuts = Dict(
        :sumE => (0, 3500),
        :e1 => (350, 5000),
        :e2 => (350, 5000),
        :dy => (0, 100),
        :dz => (0, 100),
        :trackLength1 => (400, 1000),
        :trackLength2 => (400, 1000),
        :deltaCaloTime => (0, 1.5),
        :phi => (0,180)
    )



    simu_roi_1 = (
        sumEsimu = data_cuts[:sumE],
        simulatedEnergy1 = data_cuts[:e1],
        simulatedEnergy2 = data_cuts[:e2],
        dy = data_cuts[:dy],
        dz = data_cuts[:dz],
        trackLength1 = data_cuts[:trackLength1],
        trackLength2 = data_cuts[:trackLength2],
        deltaCaloTime = data_cuts[:deltaCaloTime],
        phi = data_cuts[:phi]
    )

    data_roi_1 = (
        sumE = data_cuts[:sumE],
        reconstructedEnergy1 = data_cuts[:e1],
        reconstructedEnergy2 = data_cuts[:e2],
        dy = data_cuts[:dy],
        dz = data_cuts[:dz],
        trackLength1 = data_cuts[:trackLength1],
        trackLength2 = data_cuts[:trackLength2],
        deltaCaloTime = data_cuts[:deltaCaloTime],
        phi = data_cuts[:phi]
    )

    simu_varNames_1 = keys(simu_roi_1) |> collect .|> string
    simu_p1_processes_roi_1 = load_sim_processes_ND(p1_dir, simu_varNames_1, simu_roi_1; timeMeas = p1_duration_seconds)
    simu_p2_processes_roi_1 = load_sim_processes_ND(p1_dir, simu_varNames_1, simu_roi_1; timeMeas = p2_duration_seconds)
    simu_p3_processes_roi_1 = load_sim_processes_ND(p1_dir, simu_varNames_1, simu_roi_1; timeMeas = p3_duration_seconds)


    d1_roi_1 = filter_data!(d1, data_roi_1)
    d2_roi_1 = filter_data!(d2, data_roi_1)
    d3_roi_1 = filter_data!(d3, data_roi_1)
    

    h1_data_p1 = Hist1D(d1_roi_1.sumE; binedges = (0:100:4000))
    h1_data_p2 = Hist1D(d2_roi_1.sumE; binedges = (0:100:4000))
    h1_data_p3 = Hist1D(d3_roi_1.sumE; binedges = (0:100:4000))

    # h_data_tot = h1_data_p1 + h1_data_p2 + h1_data_p3
    h_data_tot = h1_data_p2 + h1_data_p3

    h2_simu_bb_roi_1_p1 = get_expected_simu_counts(simu_p1_processes_roi_1[1], :sumEsimu, 0.12)
    h2_simu_bb_roi_1_p2 = get_expected_simu_counts(simu_p2_processes_roi_1[1], :sumEsimu, 0.12)
    h2_simu_bb_roi_1_p3 = get_expected_simu_counts(simu_p3_processes_roi_1[1], :sumEsimu, 0.12)

    h2_simu_bb_roi_1 = h2_simu_bb_roi_1_p2 + h2_simu_bb_roi_1_p3
    # h2_simu_bb_roi_1 = h2_simu_bb_roi_1_p1 + h2_simu_bb_roi_1_p2 + h2_simu_bb_roi_1_p3

    # tot_duration = p1_duration_seconds + p2_duration_seconds + p3_duration_seconds
    tot_duration = p2_duration_seconds + p3_duration_seconds
    tot_duration_days = tot_duration / (3600*24)
    let 
        eff_bb = length(simu_p1_processes_roi_1[1].data) / simu_p1_processes_roi_1[1].nTotalSim
        f = Figure(size = (1600, 800), fontsize = 30)
        a1 = Axis(f[1,1], xlabel = "sumE", ylabel = "counts per 100keV", title = "phase 2+3 data vs simu, duration = $(round(tot_duration_days, digits=2)) days")

        p1 = stephist!(a1, h_data_tot, label = "data, n = $(round(integral(h_data_tot)))", color = :black)
        errorbars!(a1, h_data_tot, whiskerwidth = 8, color = :black)
        p2 = stephist!(a1, h2_simu_bb_roi_1, label = "2nubb simu, n = $(round(integral(h2_simu_bb_roi_1))), eff = $(round(eff_bb, digits=2))", color = :blue)

        axislegend(a1)

        el = [PolyElement(color = :transparent, strokewidth = 0)]
        cuts = join(sort(string.(keys(data_cuts)) .* " in " .* string.(values(data_cuts))), "\n")
        Legend(f[1,2], el, [cuts])
        safesave("scripts/phd_final/08_data/figs/sim_vs_data.png", f)
        f
    end

end