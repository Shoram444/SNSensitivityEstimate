import Base: *, /, sqrt

function /(x::Real, h::Hist2D)
    mat = zeros(size(bincounts(h)))

    for i in eachindex(bincounts(h))
        if (bincounts(h)[i] == 0)
            continue
        else
            mat[i] = x / bincounts(h)[i]
        end
    end
    h.bincoutns .= mat
    return h
end

function /(h::Hist2D, x::Real) 
    mat = zeros(size(bincounts(h)))

    for i in eachindex(bincounts(h))
        if (bincounts(h)[i] == 0)
            continue
        else
            mat[i] =  bincounts(h)[i] / x
        end
    end
    h.bincoutns .= mat
    return h
end


# function *(x::Real, h::Hist2D)
#     return *(h, x)
# end

function sqrt(h::Hist2D) 
    mat = zeros(size(bincounts(h)))

    for i in eachindex(bincounts(h))
        mat[i] = sqrt(bincounts(h)[i]) 
    end

    h.bincounts .= mat
    return h
end


"""
    get_max_bin(h2d::Hist2D)

Returns a Dict with the following keys:
    + :minBinEdge => minBinCenterEdge - binStepHalf ->  defines the lower edge of ROI
    + :maxBinEdge => maxBinCenterEdge + binStepHalf ->  defines the upper edge of ROI
    + :maxBinCount => maxBinCount ->  provides the maximum bin counts of the 2D histogram
"""
function get_max_bin(h2d::Hist2D)
    # check whether the histogram is square (both axes must have same binning)
    if( binedges(h2d)[1].isuniform || binedges(h2d)[1].isuniform  )
        (step(binedges(h2d)[1].uniform_edges) != step(binedges(h2d)[1].uniform_edges)) && error("bins must be the same!")
        
        # if uniform edges = bin_width is same for all bins
        halfBinWidth = step(binedges(h2d)[1].uniform_edges) / 2.
    elseif( !binedges(h2d)[1].isuniform || !binedges(h2d)[1].isuniform )
        ((binedges(h2d)[1].nonuniform_edges) != (binedges(h2d)[1].nonuniform_edges)) && error("bins must be the same!")
    
        diffs = diff(binedges(h2d)[1].nonuniform_edges)
        halfBinWidths = diffs ./ 2.
    end

    if( binedges(h2d)[1].isuniform )
        BinID = argmax(bincounts(h2d))
        minBinCenter, maxBinCenter = bincenters(h2d)[1][BinID[1]], bincenters(h2d)[1][BinID[2]]
        maxBinCount = lookup(h2d, minBinCenter, maxBinCenter)

        return Dict(
            :minBinEdge => minBinCenter - halfBinWidth,
            :maxBinEdge => maxBinCenter + halfBinWidth,
            :maxBinCount => maxBinCount
        )
    else
        BinID = argmax(bincounts(h2d))
        minBinCenter, maxBinCenter = bincenters(h2d)[1][BinID[1]], bincenters(h2d)[1][BinID[2]]
        maxBinCount = lookup(h2d, minBinCenter, maxBinCenter)
        return Dict(
            :minBinEdge => minBinCenter - halfBinWidths[BinID[1]],
            :maxBinEdge => maxBinCenter + halfBinWidths[BinID[2]],
            :maxBinCount => maxBinCount
        )
    end
end

function FHist.lookup(weights::Matrix, x::Real, y::Real, binning)
    h2d = Hist2D(Float64, bins=(binning, binning))
    h2d.hist.weights = [bc.val for bc in weights]

    rx, ry = binedges(h2d)
    !(first(rx) <= x <= last(rx)) && return missing
    !(first(ry) <= y <= last(ry)) && return missing

    return weights[FHist._edge_binindex(rx, x), FHist._edge_binindex(ry, y)]
end

FHist.binedges(h::Histogram) = h.edges[1]


function get_pseudo_spectrum(h::Hist1D)
    data_bkg = [first(FHist.sample(h)) for i=1:rand(Poisson(round(Int, integral(h))))] 
    Hist1D( Int.(data_bkg); counttype = Int, binedges= binedges(h) |> collect .|> Int )
end