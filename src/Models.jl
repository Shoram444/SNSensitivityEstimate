"""
    `turing_normal_exponential_model(d, a, b; FWHM=0.08, Q=SNparams["Q"])`
    
    # Mixutre model built with normal and exponential components
    - function parameters:
      - `d` - data; **must be shifted to 0 for proper results!** 
      - `a` - xoffset (used to shift the Q-value); equivalent to lower boundary of fitting region
      - `b` - upper bound of fitting region
      - `FWHM` - resolution in [% @ 1MeV]; used to calculate sigma at Q_keV
      - `Q` - Q-value of the normally distributed signal

    - priors:
      - `lambda`: decay constant (watch out, julia uses inverse of lambda)
      - `nSig`: number of signal events
      - `nBkg`: number of background events

"""
@model function turing_normal_exponential_model(d, a, b; FWHM=0.08, Q=SNparams["Q"])
    mu          = Q -a
    sigma_at_Q  = get_sigma_keV(Q, FWHM)
    n           = length(d)

    # Priors 
    lambda      ~ Uniform(1e-5,1e3)
    nSig        ~ Uniform(1e-15, 10)
    nBkg        ~ Uniform(1e-5, 1e3)

    # Likelihood
    l1 = nSig*pdf.(Normal(mu,sigma_at_Q), d) + nBkg*pdf.(Exponential(lambda), d) ./ (nSig+nBkg)

    l1 = sum(log.(l1))
    l2 = logpdf(Poisson(nSig+nBkg), n)

    return Turing.@addlogprob! l1+l2
end



####### 
## BAT

function BAT_likelihood_extended(NtotExpected::Real, Nobs::Real)
    pdf(Poisson(NtotExpected), Nobs)
end


"""
    BAT_normal_exponential_model

    fields: 
    - data (shifted by offset?)
    - Q
    - FWHM
"""
mutable struct BAT_normal_exponential_model
    data
    Q
    FWHM
end
@inline DensityInterface.DensityKind(::BAT_normal_exponential_model) = IsDensity()

function BAT_likelihood_normal_exponential_model(nB::Real, nS::Real, lambda::Real, x::Real; Q=2997.0, sigma = 60.0)
    (nB * pdf(Exponential(lambda), x)  +
     nS * pdf(Normal(Q, sigma),  x)) *
    inv( nB + nS ) 
end


function DensityInterface.logdensityof(m::BAT_normal_exponential_model, p)
    n = length(m.data)
    nB, nS, lambda = p

    sigma = get_sigma_keV(m.Q, m.FWHM)
    l_a = BAT_likelihood_normal_exponential_model.(nB, nS, lambda, m.data; Q=m.Q, sigma=sigma)

    ll_a = sum(log.(l_a))
    ll_b = log(BAT_likelihood_extended(nB+nS, n))
    ll_a + ll_b
end



##### Binned model
my_pdf(h::Hist1D, x::Real) = ismissing(lookup(h, x)) ? eps(Float64) : lookup(h, x)


function log_pdf_poisson(λ::T, k::U) where {T<:Real,U<:Real}
    R = float(promote_type(T,U))
    if λ >= 0 && k >= 0 && isinteger(k)
        a = iszero(k) ? R(k) : R(log(λ) * k)
        lg = logabsgamma(k + 1)[1]
        R(a - λ - lg)
    else
        R(-Inf)
    end
end

function make_hist_likelihood_dirichlet(h::Hist1D, f::Function)
    logfuncdensity(function(p)
        # Histogram counts for each bin as an array:
        counts = bincounts(h)
        n = sum(counts)
        # Histogram binning, has length (length(counts) + 1):
        binning = binedges(h)

        # sum log-likelihood value over bins:
        ll_value::Float64 = eps(Float64)
        
        for i in eachindex(counts)

            # Get information about current bin:
            bin_left, bin_right = binning[i], binning[i+1]
            bin_width = bin_right - bin_left
            bin_center = (bin_right + bin_left) / 2

            observed_counts = counts[i]

            # Simple mid-point rule integration of fit function `f` over bin:
            expected_counts = max(bin_width * f(Vector(p.w), bin_center) * n, eps(Float64))

            # expected_counts = bin_width * f(Vector(p.w), bin_center) * n
            # expected_counts == 0.0 && continue

            log_pdf_poisson(expected_counts, observed_counts) == -Inf && continue
            # Add log of Poisson probability for bin:
            ll_value += log_pdf_poisson(expected_counts, observed_counts)
            
        end
        ll_value
        return ll_value
    end )
end

function make_hist_likelihood_uniform(h::Hist1D, f::Function)
    logfuncdensity(function(p)
        # Histogram counts for each bin as an array:
        counts = bincounts(h)
        n = sum(counts)
        # Histogram binning, has length (length(counts) + 1):
        binning = binedges(h)

        # sum log-likelihood value over bins:
        ll_value::Float64 = eps(Float64)
        
        for i in eachindex(counts)

            # Get information about current bin:
            bin_left, bin_right = binning[i], binning[i+1]
            bin_width = bin_right - bin_left
            bin_center = (bin_right + bin_left) / 2

            observed_counts = counts[i]

            # Simple mid-point rule integration of fit function `f` over bin:
            expected_counts = max(bin_width * f(p, bin_center) * n, eps(Float64))

            # expected_counts = bin_width * f(p, bin_center) * n
            # expected_counts == 0.0 && continue

            log_pdf_poisson(expected_counts, observed_counts) == -Inf && continue
            # Add log of Poisson probability for bin:
            ll_value += log_pdf_poisson(expected_counts, observed_counts)
            
        end
        ll_value
        return ll_value
    end )
end

function f_dirichlet(pars::Vector{Float64}, x::Real, s_hist::Hist1D, b_hists::Vector{<:Hist1D})
    if length(pars) != length(b_hists) + 1
        error("Number of parameters must be equal to number of histograms + 1")
    end
    total_rate = sum(pars)
    sig = pars[1] / total_rate * s_hist
    bkg = [ pars[i] / total_rate * b_hists[i-1] for i in 2:length(pars)]
    th = normalize( sum(vcat(sig, bkg)), width = true )
    
    return my_pdf(th, x) 
end

function f_dirichlet(pars::Vector{Float64}, x::Real, s_hist::Hist1D, b_hist::Hist1D)
    total_rate = sum(pars)
    sig = pars[1] / total_rate * s_hist
    bkg = pars[2] / total_rate * b_hist
    th = normalize( sum(vcat(sig, bkg)), width = true )
    
    return my_pdf(th, x) 
end


function f_uniform_bkg(pars::NamedTuple{(:As, :Ab)}, x::Real, s_hist::Hist1D, b_hists::Vector{<:Hist1D})
    if length(pars.Ab) != length(b_hists)
        error("Number of parameters must be equal to number of histograms + 1")
    end
    As = pars.As
    Ab = pars.Ab |> collect
    total_rate = As + sum(Ab)

    sig = (pars.As / total_rate )*s_hist
    bkg = [ (pars.Ab[i] / total_rate)* b_hists[i] for i in 1:length(pars.Ab)]
    th = normalize( sum(vcat(sig, bkg)), width = true )
    
    return my_pdf(th, x) 
end

function f_uniform_bkg(pars::NamedTuple{(:As, :Ab)}, x::Real, s_hist::Hist1D, b_hists::Hist1D)
    Ab = if typeof(pars.Ab) == Float64
        pars.Ab
    else
        pars.Ab[1]
    end
    total_rate = pars.As + Ab
    sig = pars.As / total_rate * s_hist
    bkg = Ab / total_rate * b_hists 

    th = sum(vcat(sig, bkg))
    
    if sum(bincounts(th)) == 0
        return eps(Float64)  # Avoid zero probabilities
    end

    th = normalize(th, width = true)
    return my_pdf(th, x)
end


function get_sens_bayes_uniform(bkg_hist::Vector{<:Hist1D}, f2, signal, prior; ROI_a = ROI_a, ROI_b = ROI_b, nsteps = 10^4, nchains = 4)
    sample_hists = [get_pseudo_spectrum(b) for b in bkg_hist] 
    data_hist = merge(sample_hists...)

    my_likelihood = make_hist_likelihood_uniform(data_hist, f2)
    
    burnin = MCMCMultiCycleBurnin(max_ncycles = 50, nsteps_final=1000)
    mcmcalgo = MetropolisHastings(weighting = RepetitionWeighting(), tuning = AdaptiveMHTuning())

    posterior = PosteriorMeasure(my_likelihood, prior)
    samples, _ = bat_sample( posterior, MCMCSampling(mcalg = mcmcalgo, burnin = burnin, nsteps = nsteps, nchains = nchains))

    binned_unshaped_samples, _ = bat_transform(Vector, samples)

    nDataPoints = integral(data_hist)
    muS = [par[1] for par in binned_unshaped_samples.v]
    exp_mu_signal_90 = quantile( muS,0.9) * nDataPoints
    Na = 6.02214e23
    m = 6.067
    t = 2.88
    W = 0.08192
    eff = lookup(signal, ROI_a, ROI_b)
    t12 = log(2) * (Na * m * t * eff / W) / exp_mu_signal_90
    return t12
end