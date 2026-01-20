using DrWatson
@quickactivate("SNSensitivityEstimate")

using SNSensitivityEstimate, UnROOT, CairoMakie, LaTeXStrings, StatsBase
using Measurements
using FHist, DataFramesMeta

include(srcdir("params/Params.jl"))

begin
    duration_p0 = 38.8*24*3600 # s
    duration_p1 =  26*24*3600 # s
    tot_dur = duration_p0 + duration_p1

    f1 = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/output.root")
    d1 = LazyTree(f1, "tree", keys(f1["tree"])) |> DataFrame

    f0 = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/output.root")
    d0 = LazyTree(f0, "tree", keys(f0["tree"])) |> DataFrame

    begin
        # pint_cut = 0.02
        # pext_cut_low = 0.00 
        # pext_cut_high = 0.01

        deltaCaloTime = 1.5

        phi_low, phi_high = 10, 180 # degrees
        E_cut = 350 # keV
        dy_cut = 100 # mm
        dz_cut = 100 # mm

        l_min, l_max = 0, 5000
    end

    params = Dict(
        # "pint_cut" => 0.04,
        # "pext_cut_low" => 0.00,
        # "pext_cut_high" => 0.01,
        "deltaCaloTime" => deltaCaloTime,
        "phi_low" => phi_low,
        "phi_high" => phi_high,
        "E_cut" => E_cut,
        "dy_cut" => dy_cut,
        "dz_cut" => dz_cut,
        "l_min" => l_min,
        "l_max" => l_max
    )

    begin
        d11 = @chain d1 begin
            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)
        end

        d12 = @chain d11 begin
            @subset :reconstructedEnergy1 .> E_cut 
            @subset :reconstructedEnergy2 .> E_cut
            @subset phi_low .< :phi .< phi_high
        end

        d13 = @chain d12 begin
            @subset :dy .< dy_cut
            @subset :dz .< dz_cut
        end

        d14 = @chain d13 begin
            @subset  :deltaCaloTime .<= deltaCaloTime
            @subset l_min .< :trackLength1 .< l_max
            @subset l_min .< :trackLength2 .< l_max
        end

    end

    begin
        d01 = @chain d0 begin
            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)
        end

        d02 = @chain d01 begin
            @subset :reconstructedEnergy1 .> E_cut 
            @subset :reconstructedEnergy2 .> E_cut
            @subset phi_low .< :phi .< phi_high
        end

        d03 = @chain d02 begin
            @subset :dy .< dy_cut
            @subset :dz .< dz_cut
        end

        d04 = @chain d03 begin
            @subset  :deltaCaloTime .<= deltaCaloTime
            @subset l_min .< :trackLength1 .< l_max
            @subset l_min .< :trackLength2 .< l_max

        end

    end

    let 
        f = Figure(size = (800, 600), fontsize = 20)
        a = Axis(f[1, 1]; title = "combined phase0 + phase1",xlabel = L"E_{sum} (keV) $$", ylabel = "Counts", yscale = log10, limits = (0,3500, 1, 1e6))
        stephist!(a, Hist1D(d01.esum; binedges=0:100:5000) + Hist1D(d11.esum; binedges=0:100:5000), label = "2e topology", color = :blue, linewidth = 2)
        errorbars!(a, Hist1D(d01.esum; binedges=0:100:5000) + Hist1D(d11.esum; binedges=0:100:5000); color = :blue, whiskerwidth = 8)

        stephist!(a, Hist1D(d02.esum; binedges=0:100:5000) + Hist1D(d12.esum; binedges=0:100:5000), label = L"+ E_i > %$(E_cut) $$", color = :red, linewidth = 2)
        errorbars!(a, Hist1D(d02.esum; binedges=0:100:5000) + Hist1D(d12.esum; binedges=0:100:5000); color = :red, whiskerwidth = 8)

        stephist!(a, Hist1D(d03.esum; binedges=0:100:5000) + Hist1D(d13.esum; binedges=0:100:5000), label = L"+ |\Delta y| < %$(dy_cut) mm, |\Delta z| < %$(dz_cut) mm $$", color = :green, linewidth = 2)
        errorbars!(a, Hist1D(d03.esum; binedges=0:100:5000) + Hist1D(d13.esum; binedges=0:100:5000); color = :green, whiskerwidth = 8)

        stephist!(a, Hist1D(d04.esum; binedges=0:100:5000) + Hist1D(d14.esum; binedges=0:100:5000), label = L"+ |\Delta t_{calo}| < %$(deltaCaloTime) ns $$", color = :black, linewidth = 2)
        errorbars!(a, Hist1D(d04.esum; binedges=0:100:5000) + Hist1D(d14.esum; binedges=0:100:5000); color = :black, whiskerwidth = 8)

        axislegend(a; position = :rt)
        save("/home/maros/Work/Phd/SNSensitivityEstimate/scripts/data/combined_p0_p1/combined_p0_p1_cuts.png", f, px_per_unit = 5)
        f
        
    end

    begin
        fsim_bb_0 = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_bb_foil_bulk.root")
        dsim_bb_0 = LazyTree(fsim_bb_0, "tree", keys(fsim_bb_0["tree"])) |> DataFrame

        dsim_bb_0 = @chain dsim_bb_0 begin
            @subset :reconstructedEnergy1 .> E_cut
            @subset :reconstructedEnergy2 .> E_cut
            @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset phi_low .< :phi .< phi_high

            @subset :dy .< dy_cut
            @subset :dz .< dz_cut

            @subset  :deltaCaloTime .<= deltaCaloTime
            @subset l_min .< :trackLength1 .< l_max
            @subset l_min .< :trackLength2 .< l_max

        end
        eff_bb_0 = nrow(dsim_bb_0) / (50e3*25)
        act = halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"])
        m = 6.11
        n =  eff_bb_0 * act * m * (24*3600) * (duration_p0/tot_dur) # number of events in p0 duration scaled to 1 year

        hsim_bb_0 = normalize(Hist1D(dsim_bb_0.esum; binedges=0:100:5000), width = false) * n # normalized to 1e6 BB decays

        fsim_bb_1 = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/sim/mva_bb_foil_bulk.root")
        dsim_bb_1 = LazyTree(fsim_bb_1, "tree", keys(fsim_bb_1["tree"])) |> DataFrame

        dsim_bb_1 = @chain dsim_bb_1 begin
            @subset :reconstructedEnergy1 .> E_cut
            @subset :reconstructedEnergy2 .> E_cut
            @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset phi_low .< :phi .< phi_high

            @subset :dy .< dy_cut
            @subset :dz .< dz_cut

            @subset  :deltaCaloTime .<= deltaCaloTime
            @subset l_min .< :trackLength1 .< l_max
            @subset l_min .< :trackLength2 .< l_max

        end
        eff_bb_1 = nrow(dsim_bb_1) / (50e3*25)
        act = halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"])
        m = 6.11
        n =  eff_bb_1 * act * m * (24*3600) * (duration_p1/tot_dur) # number of events in p0 duration scaled to 1 year

        hsim_bb_1 = normalize(Hist1D(dsim_bb_1.esum; binedges=0:100:5000), width = false) * n # normalized to 1e6 BB decays
        hsim_bb = hsim_bb_0 + hsim_bb_1
    end

    begin
        d_radon_0 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_Bi214_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max
            end

        end
        eff_radon_0 = nrow(d_radon_0) / (250e3*25)
        act_radon_0 = 0.095 # Bq/m3
        n_radon_0 =  eff_radon_0 * act_radon_0 * SNparams["gasVolume"] * 24*3600 * (duration_p0/tot_dur)

        hsim_radon_0 = normalize(Hist1D(d_radon_0.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsim_radon_0 = hsim_radon_0 * (n_radon_0)

        d_radon_1 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/sim/mva_Bi214_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max
            end

        end
        eff_radon_1 = nrow(d_radon_1) / (250e3*25)
        act_radon_1 = 0.095 # Bq/m3
        n_radon_1 =  eff_radon_1 * act_radon_1 * SNparams["gasVolume"] * 24*3600 * (duration_p1/tot_dur)

        hsim_radon_1 = normalize(Hist1D(d_radon_1.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsim_radon_1 = hsim_radon_1 * (n_radon_1)

        hsim_radon = hsim_radon_0 + hsim_radon_1
    end

    begin
        d_Pa234m_0 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_Pa234m_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max

            end
        end

        eff_Pa_0 = nrow(d_Pa234m_0) / (25*50e3)
        act_Pa_0 = BkgActivityParams[:Pa234m_foil_bulk]
        n_Pa_0 =  eff_Pa_0 * act_Pa_0 * SNparams["foilMass"] * 24*3600 * (duration_p0/tot_dur)

        hsim_Pa_0 = normalize(Hist1D(d_Pa234m_0.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsim_Pa_0 = hsim_Pa_0 * (n_Pa_0)

            d_Pa234m_1 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_Pa234m_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max

            end
        end

        eff_Pa_1 = nrow(d_Pa234m_1) / (25*50e3)
        act_Pa_1 = BkgActivityParams[:Pa234m_foil_bulk]
        n_Pa_1 =  eff_Pa_1 * act_Pa_1 * SNparams["foilMass"] * 24*3600 * (duration_p1/tot_dur)

        hsim_Pa_1 = normalize(Hist1D(d_Pa234m_1.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsim_Pa_1 = hsim_Pa_1 * (n_Pa_1)
        hsim_Pa = hsim_Pa_0 + hsim_Pa_1
    end

    begin
        d_K40_0 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_K40_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max

            end
        end
        eff_K_0 = nrow(d_K40_0) / (50e3*24)
        act_K_0 = BkgActivityParams[:K40_foil_bulk]
        n_K_0 =  eff_K_0 * act_K_0 * SNparams["foilMass"] * 24*3600 * (duration_p0/tot_dur)

        hsim_K_0 = normalize(Hist1D(d_K40_0.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsim_K_0 = hsim_K_0 * (n_K_0)

        d_K40_1 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/sim/mva_K40_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max

            end
        end
        eff_K_1 = nrow(d_K40_1) / (50e3*24)
        act_K_1 = BkgActivityParams[:K40_foil_bulk]
        n_K_1 =  eff_K_1 * act_K_1 * SNparams["foilMass"] * 24*3600 * (duration_p1/tot_dur)

        hsim_K_1 = normalize(Hist1D(d_K40_1.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsim_K_1 = hsim_K_1 * (n_K_1)   

        hsim_K = hsim_K_0 + hsim_K_1
    end

    begin
        d_Bi214_0 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_Bi214_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max
            end
        end
        effBi214_0 = nrow(d_Bi214_0) / (250e3*25)
        actBi214_0 = BkgActivityParams[:Bi214_foil_bulk]
        nBi214_0 =  effBi214_0 * actBi214_0 * SNparams["foilMass"] * 24*3600 * (duration_p0/tot_dur)

        hsimBi214_0 = normalize(Hist1D(d_Bi214_0.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsimBi214_0 = hsimBi214_0 * (nBi214_0)


        d_Bi214_1 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/sim/mva_Bi214_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max
            end
        end
        effBi214_1 = nrow(d_Bi214_1) / (250e3*25)
        actBi214_1 = BkgActivityParams[:Bi214_foil_bulk]
        nBi214_1 =  effBi214_1 * actBi214_1 * SNparams["foilMass"] * 24*3600 * (duration_p1/tot_dur)

        hsimBi214_1 = normalize(Hist1D(d_Bi214_1.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsimBi214_1 = hsimBi214_1 * (nBi214_1)

        hsimBi214 = hsimBi214_0 + hsimBi214_1
    end

    begin
        d_Tl208_0 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase0_Xalbat_300keV/sim/mva_Tl208_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max
            end
        end
        effTl208_0 = nrow(d_Tl208_0) / (50e3*24)
        actTl208_0 = BkgActivityParams[:Tl208_foil_bulk]
        nTl208_0 =  effTl208_0 * actTl208_0 * SNparams["foilMass"] * 24*3600 * (duration_p0/tot_dur)

        hsimTl208_0 = normalize(Hist1D(d_Tl208_0.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsimTl208_0 = hsimTl208_0 * (nTl208_0)


        d_Tl208_1 = let 
            f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/sim/mva_Tl208_foil_bulk.root"))
            d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
            @chain d begin
                @subset :reconstructedEnergy1 .> E_cut
                @subset :reconstructedEnergy2 .> E_cut
                @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)


                @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
                @subset phi_low .< :phi .< phi_high

                @subset :dy .< dy_cut
                @subset :dz .< dz_cut

                @subset  :deltaCaloTime .<= deltaCaloTime
                @subset l_min .< :trackLength1 .< l_max
                @subset l_min .< :trackLength2 .< l_max
            end
        end
        effTl208_1 = nrow(d_Tl208_1) / (50e3*24)
        actTl208_1 = BkgActivityParams[:Tl208_foil_bulk]
        nTl208_1 =  effTl208_1 * actTl208_1 * SNparams["foilMass"] * 24*3600 * (duration_p1/tot_dur)

        hsimTl208_1 = normalize(Hist1D(d_Tl208_1.esum; binedges=0:100:5000), width = false)  # normalized to 1e6 BB decays
        hsimTl208_1 = hsimTl208_1 * (nTl208_1) 
        
        hsimTl208 = hsimTl208_0 + hsimTl208_1
    end



    h0 = Hist1D(d04.esum; binedges=0:100:5000)
    h1 = Hist1D(d14.esum; binedges=0:100:5000)

    let 
        f = Figure(size = (800, 600), fontsize = 20)
        a = Axis(f[1, 1]; xlabel = L"E_{sum} (keV)", ylabel = "Counts")
        stephist!(a, h0, label = "Phase 0", color = :blue, linewidth = 2)
        errorbars!(a, h0; color = :blue)
        stephist!(a, h1, label = "Phase 1", color = :red, linewidth = 2)
        errorbars!(a, h1; color = :red)

        axislegend(a; position = :rt)
        f

    end

    h0perday = Hist1D(d04.esum; binedges=0:100:5000) * inv(duration_p0 / 24 / 3600) 
    h1perday = Hist1D(d14.esum; binedges=0:100:5000) * inv(duration_p1 / 24 / 3600)

    # tot_dur = (duration_p0 + duration_p1)/24/3600



    hs = [hsim_bb, hsim_radon, hsim_Pa, hsim_K, hsimBi214, hsimTl208] |> reverse
    labels = [L"MC: $2\nu\beta\beta$", L"MC: radon $$", L"MC: ${}^{234m}$Pa", L"MC: ${}^{40}$K", L"MC: ${}^{214}$Bi", L"MC: ${}^{208}$Tl"] |> reverse
    colors = ["#041E42","#BE4D00","#951272","#006630","#005C8A","#FFB948"] #|> reverse
    integrals = [round(integral(h),  sigdigits=3) for h in hs] 

    let 
        tot_data = restrict(h0perday + h1perday, 700, 2900)

        f = Figure(size = (800, 600), fontsize = 24)
        a = Axis(
            f[1, 1]; 
            xlabel = L"$E_1 + E_2$ (keV)", 
            ylabel = L"$\textrm{Counts/day/100keV}$ $$",
            title = L"Phase 0 + Phase 1 data vs simulated $2\nu\beta\beta$ (%$(sum(tot_dur)) days)",
            limits = (500,4000, 0, nothing),
            # limits = (500,4000, 1.5e-1, 1.5e1),
            # yscale = log10,
            )
        
        plot!(a, sum(hs), label = L"%$(labels[end]) (%$(integrals[end]) cts/day)", color = (colors[end], 1),)
        for i in length(hs)-1:-1:1
            plot!(a, sum(hs[1:i]), label = L"%$(labels[i]) (%$(integrals[i]) cts/day)", color = colors[i], )
        end

        errorbars!(a, tot_data, label = L"$$data: Phase 0 + 1 ", color = :black, linewidth = 2, whiskerwidth = 8)
        scatter!(a, midpoints(binedges(tot_data)), bincounts(tot_data); color = :black, markersize = 6)


        # plot!(a, hsim_bb, label = L"MC: $2\nu\beta\beta$", color = (:red, 0.3),)
        # plot!(a, hsim_radon, label = L"MC: radon$$", color = (:green, 0.3),)
        # plot!(a, hsim_Pa, label = L"MC: Pa234m$$", color = (:blue, 0.3),)
        # plot!(a, hsim_K, label = L"MC: K40$$", color = (:purple, 0.3),)
        # stephist!(a, hsim_bb + hsim_radon + hsim_Pa + hsim_K, label = L"MC: total $$", color = (:red, 0.9), linewidth = 4)
        axislegend(a; position = :rt)
        save("/home/maros/Work/Phd/SNSensitivityEstimate/scripts/data/combined_p0_p1/combined_p0_p1_data_vs_sim.png", f, px_per_unit = 5)
        
        sname = savename("sim_vs_data__",params, "png")
        save("/home/maros/Work/Phd/SNSensitivityEstimate/scripts/data/combined_p0_p1/figs/" * sname, f, px_per_unit = 5)

        f

    end

    # plot(normalize(Hist1D(dsim1.phi; binedges=10:10:180)))
    # plot!(normalize(Hist1D(d04.phi; binedges=10:10:180)))
    # current_figure()
end
