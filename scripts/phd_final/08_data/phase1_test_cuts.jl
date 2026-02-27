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

function get_expected_simu_counts(p, var, fwhm; bins = (0:100:4000), timeMeas = p.timeMeas)
    smeared_energy = fwhm == 0 ? getproperty.(p.data, var) : get_smeared_energy(p, var, fwhm)
    h = normalize(Hist1D(smeared_energy; binedges = bins), width = false)
    n_exp = p.activity * timeMeas * p.amount * (length(p.data) / p.nTotalSim)
    return h .* n_exp
end

function get_chi2_ndf(data, fit, nparams)
    bincount_data = bincounts(data)
    bincount_fit = bincounts(fit)

    binerrors_data = binerrors(data)
    binerrors_fit = binerrors(fit)

    chi2 = 0.0
    ndf = 0

    for i in eachindex(bincount_data)
        O = bincount_data[i]
        E = bincount_fit[i]

        s = binerrors_data[i]^2 + binerrors_fit[i]^2

        if E > 0          # avoid division by zero
            chi2 += (O - E)^2 / s
            ndf += 1
        end
    end

    ndf = ndf - nparams
    return chi2 / ndf, chi2, ndf
end


f1 = ROOTFile(datadir("data/final_phd/official_cuts/phase1_data.root"))
d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

f2 = ROOTFile(datadir("data/final_phd/official_cuts/phase2_data.root"))
d2 = LazyTree(f2, "tree", keys(f2["tree"])) |> DataFrame

f3 = ROOTFile(datadir("data/final_phd/official_cuts/phase3_data.root"))
d3 = LazyTree(f3, "tree", keys(f3["tree"])) |> DataFrame

begin
    var_data = :sumE
    var_simu = :sumEsimu
    binning = (0:100:4000)
    data_cuts = Dict(
        :sumE => (0, 3500),
        :e1 => (250, 5000),
        :e2 => (250, 5000),
        :dy => (0, 100),
        :dz => (0, 100),
        :trackLength1 => (300, 1000),
        :trackLength2 => (300, 1000),
        :deltaCaloTime => (0, 2.0),
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
    

    h1_data_p1 = Hist1D(getproperty(d1_roi_1, var_data); binedges = binning) * inv(p1_duration_days)
    h1_data_p2 = Hist1D(getproperty(d2_roi_1, var_data); binedges = binning) * inv(p2_duration_days)
    h1_data_p3 = Hist1D(getproperty(d3_roi_1, var_data); binedges = binning) * inv(p3_duration_days)

    # h_data_tot = h1_data_p1 + h1_data_p2 + h1_data_p3
    h_data_tot = h1_data_p1


    p1_K40_process = get_process("K40_foil_bulk", simu_p1_processes_roi_1) |> first
    p1_Pa234m_process = get_process("Pa234m_foil_bulk", simu_p1_processes_roi_1) |> first
    p1_bb_process = get_process("bb_foil_bulk", simu_p1_processes_roi_1) |> first
    p1_radon_process = get_process("Bi214_wire_surface", simu_p1_processes_roi_1) |> first
    # p1_K40_PMT_process = get_process("K40_PMT_glass_bulk", simu_p1_processes_roi_1) |> first

    set_activity!(p1_radon_process, 0.0625)
    set_activity!(p1_Pa234m_process, 2.85e-3)

    simu_h_p1_K40 = get_expected_simu_counts(p1_K40_process, var_simu, 0.12, bins = binning, timeMeas = (24*3600)) * 0.8
    simu_h_p1_Pa234m = get_expected_simu_counts(p1_Pa234m_process, var_simu, 0.12, bins = binning, timeMeas = (24*3600)) * 0.8
    simu_h_p1_bb = get_expected_simu_counts(p1_bb_process, var_simu, 0.12, bins = binning, timeMeas = (24*3600)) * 0.8
    # simu_h_p1_K40_PMT = get_expected_simu_counts(p1_K40_PMT_process, var_simu, 0.12, bins = binning, timeMeas = (24*3600)) * 0.8
    simu_h_p1_radon = get_expected_simu_counts(p1_radon_process, var_simu, 0.12, bins = binning, timeMeas = (24*3600)) * 0.8




    # tot_duration = p1_duration_seconds + p2_duration_seconds + p3_duration_seconds
    let 
        eff_bb = length(simu_p1_processes_roi_1[1].data) / simu_p1_processes_roi_1[1].nTotalSim
        f = Figure(size = (1600, 1200), fontsize = 30)
        a1 = Axis(f[1,1], xlabel = "sumE", ylabel = "counts / 100keV / day", title = "phase 1 data vs simu", limits = (binning[1], binning[end], 0, nothing))

        p1 = scatter!(a1, h_data_tot, label = L"data, n = %$(round(integral(h_data_tot))) evt/day $$", color = :black, markersize = 12)
        errorbars!(a1, h_data_tot, whiskerwidth = 10, color = :black)
        p2 = stephist!(a1, simu_h_p1_K40, label = L"MC $^{40}$K, a = 58.7 mBq/kg", color = :orange,linewidth =3)
        p3 = stephist!(a1, simu_h_p1_Pa234m, label = L"MC $^{234\text{m}}$Pa, a = 17.3 mBq/kg", color = :green,linewidth =3)
        p4 = stephist!(a1, simu_h_p1_bb, label = L"MC $2\nu\beta\beta$, $T^{1/2} = 8.69 \times 10^{19}$ yr", color = :blue, linewidth =3)
        p5 = stephist!(a1, simu_h_p1_radon, label = L"MC radon, a= 95 mBq/$\textrm{m^3}$", color = :purple, linewidth =3)
        # p6 = stephist!(a1, simu_h_p1_K40_PMT, label = "K40 PMT simu, n = $(round(integral(simu_h_p1_K40_PMT)))", color = :cyan, linewidth =3)
        h_simu_total = simu_h_p1_K40 + simu_h_p1_Pa234m + simu_h_p1_bb + simu_h_p1_radon # + simu_h_p1_K40_PMT
        
        stephist!(a1, h_simu_total, label = L"MC total $$", color = :red, linewidth = 4)
        errorbars!(a1, h_simu_total, whiskerwidth = 10, color = :red)

        axislegend(a1, nbands = 2, position = :rt)

        grid = GridLayout(f[1, 2], tellheight = false)

        chi_ndf = get_chi2_ndf(h_data_tot, h_simu_total, 0)

        el = [PolyElement(color = :transparent, strokewidth = 0)]
        cuts = join(sort(string.(keys(data_cuts)) .* " in " .* string.(values(data_cuts))), "\n")
        Legend(grid[1,1], el, [cuts])

        chistring = L"$\chi^2$ / ndf = \n %$(round(chi_ndf[2], digits=2)) / %$(round(chi_ndf[3], digits=2)) = \n %$(round(chi_ndf[1], digits=2))"

        Legend(grid[2,1], [PolyElement(color = :transparent, strokewidth = 0)], [chistring])
        
        # vspan!(a1, [2600], [3100], color=:black)

        a2 = Axis(f[2,1], xlabel = "sumE", ylabel = "(data - simu)", title = "residuals", limits = (binning[1], binning[end], nothing, nothing))
        
        xs = bincenters(h_data_tot)
        ys = bincounts(h_data_tot .- h_simu_total)
        
        scatter!(a2, xs, ys, label = "residuals", color = :black, markersize = 12)
        errorbars!(a2, xs, ys, binerrors(h_data_tot .- h_simu_total),whiskerwidth = 10, color = :black)
        hlines!(a2, [0], color = :red, linestyle = :dash, linewidth = 2, label = "zero line")
        axislegend(a2)

        # vspan!(a2, [2600], [3100], color=:black)

        rowsize!(f.layout, 2, Relative(0.25))
        safesave("scripts/phd_final/08_data/figs/sim_vs_data.png", f)
        f
    end

end