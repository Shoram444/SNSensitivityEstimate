
using BAT, SpecialFunctions
using Measurements

function log_pdf_poisson(λ::Real, k::Real)
    if λ > 0 && k ≥ 0
        return k*log(λ) - λ - logabsgamma(k+1)[1]
    else
        return -Inf
    end
end
function make_hist_likelihood(data_hist::Hist1D, bkg_hists::Vector{<:Hist1D})

    data_counts = bincounts(data_hist)
    nbins = length(data_counts)
    ncomp = length(bkg_hists)

    logfuncdensity(function(p)

        μ = p.μ  # vector of expected counts per component

        ll = 0.0

        for k in 1:nbins

            λk = 0.0

            for i in 1:ncomp
                λk += μ[i] * bincounts(bkg_hists[i])[k]
            end

            λk = max(λk, eps())  # numerical protection

            ll += log_pdf_poisson(λk, data_counts[k])
        end

        return ll
    end)
end




function build_fit_histograms(bkg_hists_normed, μ)
    hists = Hist1D[]


    for i in 1:length(bkg_hists_normed)
        # FHist stores Float64 bincounts, so we scale with nominal values.
        h = Measurements.value(μ[i]) * bkg_hists_normed[i]
        push!(hists, h)
    end

    return hists
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

        if E > 5         # avoid division by zero
            chi2 += (O - E)^2 / s
            ndf += 1
        end
    end

    ndf = ndf - nparams
    return chi2 / ndf, chi2, ndf
end


function make_hist_likelihood_free(data_hist, bkg_hists)

    data_counts = bincounts(data_hist)
    bkg_counts = [bincounts(h) for h in bkg_hists]

    nbins = length(data_counts)

    logfuncdensity(function(p)

        # --- Compute μ from physics ---

        μ_K40   = p.μ_K40   
        μ_Pa    = p.μ_Pa
        μ_bb    = p.μ_bb
        μ_radon = p.μ_radon
        # μ_flat  = p.μ_flat

        μ_internal = p.μ_internal
        μ_Bi210 = p.μ_Bi210
        μ_detector = p.μ_detector
        μ_external = p.μ_external

        μ_vec = (μ_bb,μ_K40, μ_Pa, μ_internal, μ_radon,  μ_Bi210, μ_detector, μ_external)

        ll = 0.0

        for k in 1:nbins

            λk = 0.0

            for i in 1:length(bkg_hists)
                λk += μ_vec[i] * bkg_counts[i][k]
            end

            λk = max(λk, eps())
            ll += log_pdf_poisson(λk, data_counts[k])
        end

        return ll
    end)
end