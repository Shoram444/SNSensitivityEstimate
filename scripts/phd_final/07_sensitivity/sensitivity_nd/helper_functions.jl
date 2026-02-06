function sum_neutron_contributions(neutron_processes, a,b)
    fast_neutrons = filter(p -> occursin("fast", p.isotopeName), neutron_processes)
    thermal_neutrons = filter(p -> occursin("thermal", p.isotopeName), neutron_processes)   
    h_fast = sum(restrict(get_bkg_counts_1D(p), a,b) for p in fast_neutrons)
    h_thermal = sum(restrict(get_bkg_counts_1D(p), a,b) for p in thermal_neutrons)
    fast_sum = sum(measurement.(bincounts(restrict(get_bkg_counts_1D(p), a, b)), binerrors(restrict(get_bkg_counts_1D(p), a, b))) for p in fast_neutrons)
    thermal_sum = sum(measurement.(bincounts(restrict(get_bkg_counts_1D(p), a, b)), binerrors(restrict(get_bkg_counts_1D(p), a, b))) for p in thermal_neutrons)
    return fast_neutrons, thermal_neutrons, h_fast, h_thermal, fast_sum, thermal_sum
end

function get_roi_bkg_counts_hist(
    p::DataProcessND, 
    roi::NamedTuple,
    bins,
    mode
)
    data = p.data
    varNames = keys(roi)
    if roi !== nothing
        for (i,n) in enumerate(varNames)
            data = filter(x -> getproperty(x, Symbol(n)) > roi[i][1] && getproperty(x, Symbol(n)) < roi[i][2], data)
        end
    end
    p.data = data
    d = getproperty.(data, mode)

    bkg_hist = normalize(Hist1D(d;binedges = bins); width = false)

    ε = get_roi_effciencyND(p, roi).eff
    # @show p.isotopeName, ε

    bkg_hist.bincounts .= bkg_hist.bincounts .* p.amount * ε * p.activity * p.timeMeas 
    return bkg_hist
end



function get_background_counts_hists(processes, roi, bins, mode)
    h_bb = Hist1D(;binedges = bins)
    h_internal = Hist1D(;binedges = bins)
    h_radon = Hist1D(;binedges = bins)
    h_Bi210 = Hist1D(;binedges = bins)
    h_detector = Hist1D(;binedges = bins)
    h_external = Hist1D(;binedges = bins)

    bb = ["bb_foil_bulk"]
    internals = ["Bi214_foil_bulk", "K40_foil_bulk", "Pa234m_foil_bulk", "Tl208_foil_bulk"]
    radon = ["Bi214_wire_surface"]
    Bi210 = ["Bi210_wire_surface"]
    detectors = ["Bi214_PMT_glass_bulk","Bi214_PMT5_glass_bulk","K40_PMT_glass_bulk","K40_PMT5_glass_bulk","Tl208_PMT_glass_bulk","Tl208_PMT5_glass_bulk","Bi214_calo_8inch_scin_bulk","K40_calo_8inch_scin_bulk", ]
    externals = ["Tl208_hall_bulk", "Bi214_hall_bulk", "K40_hall_bulk"]

    for p in processes
        h = get_roi_bkg_counts_hist(p, roi, bins, mode)

        if p.isotopeName in internals && integral(h) >0
            h_internal = h_internal + h
        elseif p.isotopeName in bb && integral(h) >0
            h_bb = h_bb + h
        elseif p.isotopeName in radon && integral(h) >0
            h_radon = h_radon + h
        elseif p.isotopeName in Bi210 && integral(h) >0
            h_Bi210 = h_Bi210 + h
        elseif p.isotopeName in detectors && integral(h) >0
            h_detector = h_detector + h
        elseif integral(h) >0
            h_external = h_external + h
        end
    end
    return h_bb, h_internal, h_radon, h_Bi210, h_detector, h_external

end

function stacked_hist!(ax, hists; labels, colors, fill_area)
    if fill_area
        hist!(ax, sum(hists),
            label = isnothing(labels) ? nothing : labels[1],
            color = colors[1],
            strokewidth = 2,
            strokecolor = :black)
        errorbars!(ax, sum(hists); whiskerwidth=7, color = :black, clamp_bincounts=true)
        
    else
        stephist!(ax, sum(hists),
            label = isnothing(labels) ? nothing : labels[1],
            color = colors[1],
            linewidth = 4)        
        errorbars!(ax, sum(hists); whiskerwidth=7, color = colors[1], clamp_bincounts=true)

    end

    for i in 2:length(hists)
        if fill_area
            hist!(ax, sum(hists[i:end]),
                label = isnothing(labels) ? nothing : labels[i],
                color = colors[i],
                strokewidth = 2,
                strokecolor = :black)
            nonzero_hists = sum(filter(h -> maximum(bincounts(h)) > 0, hists[i:end])) # when histogram is empty gives NaN errors
            errorbars!(ax, nonzero_hists; whiskerwidth=7, color = :black, clamp_bincounts=true)

        else
            stephist!(ax, sum(hists[i:end]),
                label = isnothing(labels) ? nothing : labels[i],
                color = colors[i],
                linewidth = 4)
            nonzero_hists = sum(filter(h -> maximum(bincounts(h)) > 0, hists[i:end])) # when histogram is empty gives NaN errors
            errorbars!(ax, nonzero_hists; whiskerwidth=7, color = colors[i], clamp_bincounts=true)
        end

    end
end

function find_max_bincounts(hists::Hist1D)
    if typeof(hists) <: FHist.Hist1D
        return maximum(bincounts(hists))
    end
end
