import Base: *, /

function /(x::Real, h::Hist2D)
    mat = zeros(size(bincounts(h)))

    for i in eachindex(bincounts(h))
        if (bincounts(h)[i] == 0)
            continue
        else
            mat[i] = x / bincounts(h)[i]
        end
    end
    h.hist.weights = mat
    return h
end

function /(h::Hist2D, x::Real)
    mat = zeros(size(bincounts(h)))

    for i in eachindex(bincounts(h))
        mat[i] = bincounts(h)[i] / x
    end

    h.hist.weights = mat
    return h
end

function *(x::Real, h::Hist2D)
    return *(h, x)
end


"""
    get_max_bin(h2d::Hist2D)

Returns a Dict with the following keys:
    + :minBinEdge => minBinCenterEdge - binStepHalf ->  defines the lower edge of ROI
    + :maxBinEdge => maxBinCenterEdge + binStepHalf ->  defines the upper edge of ROI
    + :maxBinCount => maxBinCount ->  provides the maximum bin counts of the 2D histogram
"""
function get_max_bin(h2d::Hist2D)
    (step(binedges(h2d)[1]) != step(binedges(h2d)[1])) && error("bins must be the same!")
    binStepHalf = step(binedges(h2d)[1]) / 2
    BinID = argmax(bincounts(h2d))
    minBinEdge, maxBinEdge = bincenters(h2d)[1][BinID[1]], bincenters(h2d)[1][BinID[2]]
    maxBinCount = lookup(h2d, minBinEdge, maxBinEdge)

    return Dict(
        :minBinEdge => minBinEdge - binStepHalf,
        :maxBinEdge => maxBinEdge + binStepHalf,
        :maxBinCount => maxBinCount
    )
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