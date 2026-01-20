using DrWatson
@quickactivate("SNSensitivityEstimate")

using SNSensitivityEstimate, UnROOT, CairoMakie, LaTeXStrings, StatsBase
using Measurements
using FHist, DataFramesMeta

include(srcdir("params/SNparams.jl"))


f = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v3_phase1_Xalbat_50keV/output.root")
d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
first(select(d, [:reconstructedEnergy1, :reconstructedEnergy2]), 10)

global new_c = 299.792458 # speed of light in mm/ns
begin
    
    pint_cut = 0.04
    pext_cut_low = 0.00 
    pext_cut_high = 0.01

    deltaCaloTime = 2

    phi_low, phi_high = 0, 180 # degrees
    E_cut = 350 # keV
    dy_cut = 30 # mm
    dz_cut = 80 # mm

    d1 = @chain d begin
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)
    end

    d2 = @chain d1 begin
        @subset :reconstructedEnergy1 .> E_cut 
        @subset :reconstructedEnergy2 .> E_cut
        # @subset :esum .> E_cut
        @subset phi_low .< :phi .< phi_high
    end

    d3 = @chain d2 begin
        @subset :dy .< dy_cut
        @subset :dz .< dz_cut
    end


    d4 = @chain d3 begin
        # @subset :Pint .>= pint_cut # t1 and t2 are in ns
        # @subset :tint .<= 10
        # @subset 400 .< :trackLength2 .< 1000 
        # @subset 400 .< :trackLength1 .< 1000
    end

    d5 = @chain d4 begin
        # @subset :Pint .>= :Pext

        # @subset pext_cut_low .<= :Pext .<= pext_cut_high # external theoretical time diff in ns
        # @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)
        @subset  :deltaCaloTime .<= deltaCaloTime
    end

    ### signal sim


    d1_sim = let
        f_sim = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v2_phase1_Xalbat_50keV/sim_bb_foil_bulk_mva.root")
        d_sim = LazyTree(f_sim, "tree", keys(f_sim["tree"])) |> DataFrame

        @chain d_sim begin
            @subset :reconstructedEnergy1 .> E_cut
            @subset :reconstructedEnergy2 .> E_cut
            @rtransform :deltaCaloTime = abs(:caloTime1 - :caloTime2)

            # @subset :phi .> 15
            # @rtransform :beta1 = sqrt(:reconstructedEnergy1*(:reconstructedEnergy1 + 2*511)) / (:reconstructedEnergy1 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            # @rtransform :beta2 = sqrt(:reconstructedEnergy2*(:reconstructedEnergy2 + 2*511)) / (:reconstructedEnergy2 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            
            # @rtransform :t1 = :trackLength1 / (:beta1*new_c) #t_1 = track_length_1/(beta_1*new_c)
            # @rtransform :t2 = :trackLength2 / (:beta2*new_c) 
            
            # @rtransform :tint = ((:time1-:t1) - (:time2-:t2))^2 #((:t1 - :t2) - (:time1 - :time2))^2 #internal_theoretical_time_diff = (t_1-t_2) - diff_time_elec;
            # @rtransform :text = ((:time2-:time1)-(:t1+:t2))^2  #(abs(:time1 - :time2) - (:t1 + :t2))^2 #external_theoretical_time_diff = abs(diff_time_elec) - (t_1+t_2) ;

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            # @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
            # @subset :reconstructedEnergy2 .> E_cut
            @subset phi_low .< :phi .< phi_high


            # @rtransform :dy = abs(:y1Reconstructed - :y2Reconstructed)
            # @rtransform :dz = abs(:z1Reconstructed - :z2Reconstructed)
            @subset :dy .< dy_cut
            @subset :dz .< dz_cut
            # @subset :Pint .>= pint_cut # t1 and t2 are in ns
            # @subset pext_cut_low .<= :Pext .<= pext_cut_high
            @subset  :deltaCaloTime .<= deltaCaloTime

        end
    end

    ### background processes

    d_K40 = let 
        f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/K40_foil_bulk.root"))
        d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
        @chain d begin
            @subset :reconstructedEnergy1 .> 0
            @subset :reconstructedEnergy2 .> 0
            # @subset :phi .> 15
            # @rtransform :beta1 = sqrt(:reconstructedEnergy1*(:reconstructedEnergy1 + 2*511)) / (:reconstructedEnergy1 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            # @rtransform :beta2 = sqrt(:reconstructedEnergy2*(:reconstructedEnergy2 + 2*511)) / (:reconstructedEnergy2 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            
            # @rtransform :t1 = :trackLength1 / (:beta1*new_c) #t_1 = track_length_1/(beta_1*new_c)
            # @rtransform :t2 = :trackLength2 / (:beta2*new_c) 
            
            # @rtransform :tint = ((:time1-:t1) - (:time2-:t2))^2 #((:t1 - :t2) - (:time1 - :time2))^2 #internal_theoretical_time_diff = (t_1-t_2) - diff_time_elec;
            # @rtransform :text = ((:time2-:time1)-(:t1+:t2))^2  #(abs(:time1 - :time2) - (:t1 + :t2))^2 #external_theoretical_time_diff = abs(diff_time_elec) - (t_1+t_2) ;

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
            # @subset :reconstructedEnergy2 .> E_cut
            @subset phi_low .< :phi .< phi_high

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset :esum .> E_cut

            # @rtransform :dy = abs(:y1Reconstructed - :y2Reconstructed)
            # @rtransform :dz = abs(:z1Reconstructed - :z2Reconstructed)
            @subset :dy .< dy_cut
            @subset :dz .< dz_cut
            @subset :Pint .>= pint_cut # t1 and t2 are in ns
            @subset pext_cut_low .<= :Pext .<= pext_cut_high
        end
    end

    d_Pa234m = let 
        f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/Pa234m_foil_bulk.root"))
        d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
        @chain d begin
            @subset :reconstructedEnergy1 .> E_cut
            @subset :reconstructedEnergy2 .> E_cut
            # @subset :phi .> 15
            # @rtransform :beta1 = sqrt(:reconstructedEnergy1*(:reconstructedEnergy1 + 2*511)) / (:reconstructedEnergy1 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            # @rtransform :beta2 = sqrt(:reconstructedEnergy2*(:reconstructedEnergy2 + 2*511)) / (:reconstructedEnergy2 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            
            # @rtransform :t1 = :trackLength1 / (:beta1*new_c) #t_1 = track_length_1/(beta_1*new_c)
            # @rtransform :t2 = :trackLength2 / (:beta2*new_c) 
            
            # @rtransform :tint = ((:time1-:t1) - (:time2-:t2))^2 #((:t1 - :t2) - (:time1 - :time2))^2 #internal_theoretical_time_diff = (t_1-t_2) - diff_time_elec;
            # @rtransform :text = ((:time2-:time1)-(:t1+:t2))^2  #(abs(:time1 - :time2) - (:t1 + :t2))^2 #external_theoretical_time_diff = abs(diff_time_elec) - (t_1+t_2) ;

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
            # @subset :reconstructedEnergy2 .> E_cut
            @subset phi_low .< :phi .< phi_high

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset :esum .> E_cut

            # @rtransform :dy = abs(:y1Reconstructed - :y2Reconstructed)
            # @rtransform :dz = abs(:z1Reconstructed - :z2Reconstructed)
            @subset :dy .< dy_cut
            @subset :dz .< dz_cut
            @subset :Pint .>= pint_cut # t1 and t2 are in ns
            @subset pext_cut_low .<= :Pext .<= pext_cut_high
        end
    end

    d_radon = let 
        f = ROOTFile(datadir("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v0_phase1_50keV_old_calibration/mva/Bi214_wire_surface.root"))
        d = LazyTree(f, "tree", keys(f["tree"])) |> DataFrame
        @chain d begin
            @subset :reconstructedEnergy1 .> E_cut
            @subset :reconstructedEnergy2 .> E_cut
            # @subset :phi .> 15
            # @rtransform :beta1 = sqrt(:reconstructedEnergy1*(:reconstructedEnergy1 + 2*511)) / (:reconstructedEnergy1 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            # @rtransform :beta2 = sqrt(:reconstructedEnergy2*(:reconstructedEnergy2 + 2*511)) / (:reconstructedEnergy2 + 511) # sqrt(E_1*(E_1+2*0.511))/(E_1+0.511)
            
            # @rtransform :t1 = :trackLength1 / (:beta1*new_c) #t_1 = track_length_1/(beta_1*new_c)
            # @rtransform :t2 = :trackLength2 / (:beta2*new_c) 
            
            # @rtransform :tint = ((:time1-:t1) - (:time2-:t2))^2 #((:t1 - :t2) - (:time1 - :time2))^2 #internal_theoretical_time_diff = (t_1-t_2) - diff_time_elec;
            # @rtransform :text = ((:time2-:time1)-(:t1+:t2))^2  #(abs(:time1 - :time2) - (:t1 + :t2))^2 #external_theoretical_time_diff = abs(diff_time_elec) - (t_1+t_2) ;

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
            # @subset :reconstructedEnergy2 .> E_cut
            @subset phi_low .< :phi .< phi_high

            @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
            @subset :esum .> E_cut

            # @rtransform :dy = abs(:y1Reconstructed - :y2Reconstructed)
            # @rtransform :dz = abs(:z1Reconstructed - :z2Reconstructed)
            @subset :dy .< dy_cut
            @subset :dz .< dz_cut
            @subset :Pint .>= pint_cut # t1 and t2 are in ns
            @subset pext_cut_low .<= :Pext .<= pext_cut_high
        end
    end


    begin
        eff = nrow(d1_sim) / 1e6
        act = halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"])
        m = 6.11
        d_p1  = 2.244555e6
        n = eff * act * d_p1 * m # number of events in 1 year


        fit_rates = [0.7194023684148142, 0.0005319732686368019, 0.003246943932229697, 0.2753713801588674]
        std_rates = [0.023566868040843413, 0.000771064568315688, 0.0046650013773221, 0.022790061385081533]
        
        fit_n_std = [round(r*nrow(d5), digits = 1) for r in std_rates]
        fit_n = [round(r*nrow(d5), digits = 1) for r in fit_rates]
        fit_n_w_std = [measurement(fit_n[i], fit_n_std[i]) for i in 1:length(fit_n)]

        f = Figure(size = (1400,1000), fontsize = 20)
        h1 = Hist1D(d1.esum; binedges = 0:100:3500)
        h2 = Hist1D(d2.esum; binedges = 0:100:3500)
        h3 = Hist1D(d3.esum; binedges = 0:100:3500)
        h4 = Hist1D(d4.esum; binedges = 0:100:3500)
        h5 = Hist1D(d5.esum; binedges = 0:100:3500)


        a= Axis(f[1,1], xlabel = "E1 + E2", ylabel = "counts", title = "data bb phase1", limits = (0, 3500, nothing, nothing))

        p1 = errorbars!(a, h1, color=:black, whiskerwidth=5,)
        stairs!(a, h1, color=:black, linewidth = 2, label = L"no cuts $$")

        p2 =errorbars!(a, h2, color=Makie.wong_colors()[1], whiskerwidth=5,)
        stairs!(a, h2, color=Makie.wong_colors()[1], linewidth = 2, label = L"$\varepsilon_1: E_i >$ %$(E_cut) keV + $ \varphi \in (%$(phi_low), %$(phi_high)) {}^\circ$")

        p3 = errorbars!(a, h3, color=Makie.wong_colors()[2], whiskerwidth=5,)
        stairs!(a, h3, color=Makie.wong_colors()[2], linewidth = 2, label = L"$\varepsilon_1 + \varepsilon_2: \Delta y < $(dy_cut) mm, \Delta z < $(dz_cut) mm$")

        # p4 = errorbars!(a, h4, color=Makie.wong_colors()[3], whiskerwidth=5,)
        # stairs!(a, h4, color=Makie.wong_colors()[3], linewidth = 2, label = L"$\varepsilon_1 + \varepsilon_2 + \varepsilon_3: \Delta t_{int} <= $(pint_cut) ns^2$")

        p5= errorbars!(a, h5, color=Makie.wong_colors()[4], whiskerwidth=5,)
        stairs!(a, h5, color=Makie.wong_colors()[4], linewidth = 4, label = L"$\varepsilon_1 + \varepsilon_2 + \varepsilon_3 + \varepsilon_4: t_{ext} \in (%$(pext_cut_low), %$(pext_cut_high)) ns^2$")

        h6 = Hist1D(d1_sim.esum; binedges = 0:100:3500)
        h6 = normalize(h6, width = :false) * n
        # h6 = Hist1D(d1_sim.esum; binedges = 0:100:3500, weights = fill((nrow(d5)/nrow(d1_sim))*fit_rates[1], nrow(d1_sim)))
        # h7 = Hist1D(d_K40.esum; binedges = 0:100:3500, weights = fill((nrow(d5)/nrow(d_K40))*fit_rates[2], nrow(d_K40)))
        # h8 = Hist1D(d_Pa234m.esum; binedges = 0:100:3500, weights = fill((nrow(d5)/nrow(d_Pa234m))*fit_rates[3], nrow(d_Pa234m)))
        # h9 = Hist1D(d_radon.esum; binedges = 0:100:3500, weights = fill((nrow(d5)/nrow(d_radon))*fit_rates[4], nrow(d_radon)))

        p6 = stairs!(a, h6, color = (Makie.wong_colors()[2],1 ), label = L"$2\nu\beta\beta$ simulation", linestyle = :dash, linewidth = 3)
        # p7 = stairs!(a, h7, color = (Makie.wong_colors()[3],1 ), label = L"${}^{40}$K simulation", linestyle = :dash, linewidth = 3)
        # p8 = stairs!(a, h8, color = (Makie.wong_colors()[4],1 ), label = L"${}^{234}$Pa simulation", linestyle = :dash, linewidth = 3)
        # p10 = stairs!(a, h9, color = (Makie.wong_colors()[6],1 ), label = L"$$ radon simulation", linestyle = :dash, linewidth = 3)

        # p9 = plot!(a, h6+h7+h8+h9, color = (Makie.wong_colors()[1], 0.05), )


        # axislegend(a, position = :rt, nbanks = 2)
        Legend(
            f[1,2], 
            [
                p1,
                p2,
                p3,
                # p4,
                p5,
            ], 
            [
                L"$\varepsilon_0:$ 2$e^-$ topology only",
                L"$\varepsilon_0 + \varepsilon_1: E_i > %$(E_cut)$ keV + $\varphi \in (%$(phi_low), %$(phi_high)) {}^\circ$",
                L"$\varepsilon_0 + \varepsilon_1 + \varepsilon_2: \Delta y < %$(dy_cut) mm,~ \Delta z < %$(dz_cut) mm$",
                # L"$\varepsilon_0 + \varepsilon_1 + \varepsilon_2 + \varepsilon_3: d \in (400, 1000) mm$",
                # L"$\varepsilon_0 + \varepsilon_1 + \varepsilon_2 + \varepsilon_3: \Delta P_{int} >= %$(pint_cut)$",
                # L"$\varepsilon_0 + \varepsilon_1 + \varepsilon_2 + \varepsilon_3 + \varepsilon_4: P_{ext} \in (%$(pext_cut_low), %$(pext_cut_high))$",
                L"$\varepsilon_0 + \varepsilon_1 + \varepsilon_2 + \varepsilon_3 + \varepsilon_4: \Delta caloTime <= %$(deltaCaloTime) ns$",
            ],
            nbanks = 1,
            tellheight = false,
            tellwidth = true,
            patchsize = (40, 30)
        )

        a2 = Axis(f[2,1], xlabel = "E1 + E2", ylabel = "counts", title = "golden events: $(nrow(d5))", limits = (0, 3500, nothing, nothing))

        # plot!(a2, p5)
        # plot!(a2, p9)
        plot!(a2, p6)
        # plot!(a2, p7)
        # plot!(a2, p8)
        # plot!(a2, p10)

        errorbars!(a2, h5, color=Makie.wong_colors()[4], whiskerwidth=5,)
        stairs!(a2, h5, color=Makie.wong_colors()[4], linewidth = 4, label = L"$\varepsilon_1 + \varepsilon_2 + \varepsilon_3 + \varepsilon_4: t_{ext} \in (%$(pext_cut_low), %$(pext_cut_high)) ns^2$")


        # stairs!(a2, h5, color=Makie.wong_colors()[4], linewidth = 2, label = L"$\varepsilon_1 + \varepsilon_2 + \varepsilon_3 + \varepsilon_4: \Delta t_{ext} \in (5, 30) ns^2$")
        hidexdecorations!(a2)
        a3 = Axis(f[3,1], xlabel = "E1 + E2", ylabel = "data/MC (%)", limits = (0, 3500, 0, 200))
        # toth = h6 + h7 + h8 + h9
        toth = h6 
        h_data_bc_w_errors = measurement.(bincounts(h5), binerrors(h5))
        h_sim_bc_w_errors =  measurement.(bincounts(toth), binerrors(toth))

        ratio = h_data_bc_w_errors ./ h_sim_bc_w_errors .* 100
        
        CairoMakie.scatter!(a3, midpoints(binedges(toth)), ratio, color = :black, markersize = 10, )
        CairoMakie.errorbars!(a3, midpoints(binedges(toth)), [r.val for r in ratio],[r.err for r in ratio], color = :black, whiskerwidth = 5, label = "(data-MC)/MC (%)")
        # stairs!(a3, midpoints(binedges(h5)), bincounts(h5) ./ bincounts(h6))
        rowsize!(f.layout, 3, Relative(0.15))

        Legend(
            f[2,2], 
            [
                p5,
                p6,
                # p7,
                # p8,
                # p10,
                # p9
                
            ], 
            [
                L"data: $\varepsilon_0 + \varepsilon_1 + \varepsilon_2 + \varepsilon_3 + \varepsilon_4$",
                L"$2\nu\beta\beta$; n = %$(round(integral(h6)))",
                # L"$2\nu\beta\beta$; n = %$(fit_n_w_std[1])",
                # L"${}^{40}$K; n = %$(fit_n_w_std[2])",
                # L"${}^{234}$Pa; n = %$(fit_n_w_std[3])",
                # L"$$ radon; n = %$(fit_n_w_std[4])",
                # L"$2\nu\beta\beta + {}^{40}$K + ${}^{234}$Pa + radon simulation"
            ],
            nbanks = 1,
            tellheight = false,
            tellwidth = true,
            patchsize = (40, 30)

        )
        safesave(joinpath(scriptsdir("data"),"phase1", "sum_spectrum.png"), f, px_per_unit = 5)
        f

        # hhh = Hist1D(d5.phi; binedges = 0:5:180)
        # plot(hhh)
    end
end

let
    eff = nrow(d1_sim) / 1e6
    act = halfLife_to_activity(SNparams["Nₐ"], SNparams["W"], SNparams["SeThalf2nu"])
    m = 6.11
    d_p1  = 2.244555e6
    n = eff * act * d_p1 * m # number of events in 1 year
    f = Figure()
    a= Axis(
        f[1,1], 
        xlabel = "E1 + E2", 
        ylabel = "counts", 
        title = L"data bb phase1; \\ aquisiton time $\approx$ %$(round(d_p1/3600/24, digits =1)) days", 
        limits = (0, 3500, nothing, nothing)
    )
    stairs!(a, h5, label = "data\nn = $(Int(integral(h5)))", linewidth = 4)

    h6 = Hist1D(d1_sim.esum; binedges = 0:100:3500)
    h6 = normalize(h6, width = :false)
    plot!(a, h6 * n, label = "2nu; \nn = $(Int(round(n)))", color = (Makie.wong_colors()[2], 0.5),)
    axislegend(a, position = :rt)
    save(joinpath(scriptsdir("data/phase1"), "phase1_2e_spectrum_2nu.png"), f, px_per_unit = 5)
    f 
end


h2d_l_t = Hist2D((
    vcat(d.trackLength1, d.trackLength2),
    vcat(d.t1, d.t1)
    ); 
    binedges = (0:10:1000, 0:0.1:10),

)

plot(
    h2d_l_t,
    axis = (
            xlabel = "track_length",
            ylabel = "time"
        )
)







ff = ROOTFile("/home/maros/Work/Phd/SNSensitivityEstimate/data/data/v2_phase1_Xalbat_50keV/sim_bb_foil_bulk_mva.root")
dd = LazyTree(ff, "tree", keys(ff["tree"])) |> DataFrame

let 
    ddd = @chain dd begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset :esum .> E_cut

        @subset phi_low .< :phi .< phi_high

        # @rtransform :dy = abs(:y1Reconstructed - :y2Reconstructed)
        # @rtransform :dz = abs(:z1Reconstructed - :z2Reconstructed)

        @subset :dy .< dy_cut
        @subset :dz .< dz_cut

        # @subset 400 .< :trackLength2 .< 1000 
        # @subset 400 .< :trackLength1 .< 1000

    end
    tim = abs.(ddd.caloTime1 .- ddd.caloTime2)
    dist = ((ddd.trackLength1 ./ ddd.beta1) .+ (ddd.trackLength2 ./ ddd.beta2))

    h2d = Hist2D(
        (tim, dist);
        binedges = (0:0.1:10, 800:10:2000)
        )
    f,a,p = plot(h2d, axis = (limits =(0,10, 800, 2000),xlabel = "abs(t1-t2)", ylabel = "(l1/b1) + (l2/b2)"))
    current_figure()
    y1(x) = 228*x + 85
    y2(x) = 228*x -215
    lines!(a, 0:0.1:10, y1.(0:0.1:10), linewidth = 4)
    lines!(a, 0:0.1:10, y2.(0:0.1:10), linewidth = 4)
    f
end



let 
    ddd = @chain dd begin
        @subset :reconstructedEnergy1 .> 0
        @subset :reconstructedEnergy2 .> 0
        @rtransform :esum = :reconstructedEnergy1 + :reconstructedEnergy2
        @subset :reconstructedEnergy1 .> E_cut .|| :reconstructedEnergy2 .> E_cut
        @subset :esum .> E_cut

        @subset phi_low .< :phi .< phi_high

        # @rtransform :dy = abs(:y1Reconstructed - :y2Reconstructed)
        # @rtransform :dz = abs(:z1Reconstructed - :z2Reconstructed)

        @subset :dy .< dy_cut
        @subset :dz .< dz_cut

        # @subset 400 .< :trackLength2 .< 1000 
        # @subset 400 .< :trackLength1 .< 1000

    end
    tim = abs.((ddd.caloTime1 .- ddd.caloTime2))
    dist = abs.(((ddd.trackLength1 ./ ddd.beta1) .- (ddd.trackLength2 ./ ddd.beta2)))

    @show extrema(tim)
    @show extrema(dist)
    h2d = Hist2D(
        (tim, dist);
        binedges = (0:0.1:10, 0:5:300)
        )
    f,a,p = plot(h2d, axis = (xlabel = "abs(t1-t2)", ylabel = "abs(l1/b1 - l2/b2)"))
    current_figure()
    # y1(x) = 228*x + 85
    # y2(x) = 228*x -215
    # lines!(a, 0:0.1:10, y1.(0:0.1:10), linewidth = 4)
    # lines!(a, 0:0.1:10, y2.(0:0.1:10), linewidth = 4)
    f
end