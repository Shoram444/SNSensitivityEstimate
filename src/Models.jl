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