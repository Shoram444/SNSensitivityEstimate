@userplot StackedHist

@recipe function f(h::StackedHist; nTotSim = 1e8)
    histos = h.args[1] # container of the passed histograms
    wgts = h.args[2] # container of the passed weights

    for i in 1:length(histos)
        binedges = FHist.binedges(histos[1])
        bincounts = zeros(Float64, length(FHist.bincounts(histos[1])))

        @show length(binedges), length(bincounts)

        for j=1:i
            eff = FHist.integral(histos[j]) / nTotSim
            totalArea= FHist.integral(histos[j]) * step(binedges)
            bincounts .+= FHist.bincounts(histos[j]) .* wgts[j] ./ totalArea .* eff # Normalize to number of events and activity
        end

        @series begin
            seriestype := :stepbins
            binedges, bincounts
        end
    end
end