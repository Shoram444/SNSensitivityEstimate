Plots.@userplot StackedHist

Plots.@recipe function f(h::StackedHist)
    histos = h.args[1] # container of the passed histograms
    isylog = pop!(plotattributes, :yscale, :identity) ∈ (:log10, :log)
    the_ylims = pop!(plotattributes, :ylims, (-Inf, Inf))
    
    x = midpoints(FHist.binedges(histos[1])) 

    y = zeros(length(x), length(histos))
    fr = zeros(length(x), length(histos))

    for i in 1:length(histos)
        y[:, i] = i == 1 ? bincounts(histos[i]) : y[:, i-1] .+ bincounts(histos[i])
        fr[:, i] = i == 1 ? zeros(length(x)) : y[:, i-1]
    end

    fr[:,1] .= if isylog
        if isfinite(the_ylims[1])
            min(minimum(y) / 100, the_ylims[1])
        else
            minimum(y) / 100
        end
    else
        0
    end

    linewidth --> 0.0
    bar_width --> step(binedges(histos[1]).uniform_edges)
    fillrange  := fr
    seriestype := :bar

    x,y
end

Plots.@userplot StackedStepHist


Plots.@recipe function f(h::StackedStepHist)
    histos = h.args[1] # container of the passed histograms
    isylog = pop!(plotattributes, :yscale, :identity) ∈ (:log10, :log)
    the_ylims = pop!(plotattributes, :ylims, (-Inf, Inf))
    
    x = midpoints(FHist.binedges(histos[1])) 

    y = zeros(length(x), length(histos))    

    for i in 1:length(histos)
        y[:, i] = i == 1 ? FHist.bincounts(histos[i]) : y[:, i-1] .+ FHist.bincounts(histos[i])
    end

    seriestype := :stepmid

    x,y
end