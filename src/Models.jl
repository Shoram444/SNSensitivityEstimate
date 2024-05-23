@model function exponential_normal_model(data, a, b; FWHM=0.08, Q=SNparams["Q"])
    # first data must be offset to start at 0
    sigma_at_Q = get_sigma_keV(Q, FWHM)

    # Priors
    lambda ~ Uniform(1e-3, 1e3)
    nSignal ~ Uniform(0.0, 10.0)
    nBkg ~ Uniform(0.0, 1e3)

    # Likelihood 
    for x in data
        likelihood_normal = logpdf(Normal(Q-a, sigma_at_Q), x - a)
        likelihood_exponential = logpdf(Exponential(inv(lambda)), x - a)
        likelihood_mixture = (nSignal * likelihood_normal + nBkg * likelihood_exponential) / (nSignal + nBkg)
        likelihood_extended = logpdf(Poisson(nSignal + nBkg), length(data))
        Turing.@addlogprob! ((likelihood_extended) + (likelihood_mixture))
    end
end

@model function exponential_uniform_normal_model(data, a, b; FWHM=0.08, Q=SNparams["Q"])
    # first data must be offset to start at 0
    sigma_at_Q = get_sigma_keV(Q, FWHM)
    Q -= a

    # Priors
    lambda ~ Uniform(1e-3, 1e3)
    sigma ~ truncated(Normal(sigma_at_Q, 1), 1e-3, Inf)
    nSignal ~ Uniform(0.0, 10.0)
    nBkg1 ~ Uniform(0.0, 1e3)
    nBkg2 ~ Uniform(0.0, 1e3)

    # Likelihood 
    for x in data
        likelihood_normal = pdf(Normal(Q, sigma), x - a)
        likelihood_exponential = pdf(Exponential(inv(lambda)), x - a)
        likelihood_uniform = inv(abs(b - a))
        likelihood_mixture =
            (nSignal * likelihood_normal + nBkg1 * likelihood_exponential + nBkg2 * likelihood_uniform) *
            inv(nSignal + nBkg1 + nBkg2)
        likelihood_extended = pdf(Poisson(nSignal + nBkg1 + nBkg2), length(data))
        Turing.@addlogprob! (log(likelihood_extended) + log(likelihood_mixture))
    end
end


@model mixture_model(data, a, b; FWHM=0.08, Q=SNparams["Q"]) = begin
    # Priors
    lambda ~ Uniform(1e-3, 1e3)
    nSig ~ Uniform(0, 10)
    nBkg ~ Uniform(0, 3*length(data))
    
    # Derived parameters
    total_n = nSig + nBkg
    
    # Mixture probabilities
    p_signal = nSig / total_n
    p_background = nBkg / total_n
    
    # Likelihood
    data ~ MixtureModel([Normal(Q-a, 60), Exponential(lambda)], [p_signal, p_background])
end

"""
    `exponential_normal_mixture(d, a, b; FWHM=0.08, Q=SNparams["Q"])`
    
    # Mixutre model built with normal and exponential components
    - function parameters:
      - `d` - data; **must be shifted to 0 for proper results!** 
      - `a` - xoffset (used to shift the Q-value); equivalent to lower boundary of fitting region
      - `b` - upper bound of fitting region
      - `FWHM` - resolution in [% @ 1MeV]; used to calculate sigma at Q_keV
      - `Q` - Q-value of the normally distributed signal

    - priors:
      - `mu`: equivalent to Q value (shifted by offset a: Q-a)
      - `sigma`: energy resolution
      - `lambda`: decay constant (watch out, julia uses inverse of lambda)
      - `w`: weights equivalent to the proportion of nSig and nBkg in MixtureModel; To calculate `nSig = w[1]*lenght(data)`

"""
@model function exponential_normal_mixture(d, a, b; FWHM=0.08, Q=SNparams["Q"])
    # Priors for distributions
    mu          ~ Normal(Q - a,0.1)
    sigma_at_Q  = get_sigma_keV(Q, FWHM)
    sigma       ~ LogNormal(log(sigma_at_Q) - 0.5)
    lambda      ~ Uniform(1e-3, 1e3)

    # Priors for nSig, nBkg  
    K = 2 # number of distributions in likelihood
    w ~ Dirichlet(K, 1.0)

    # Likelihood
    dist = MixtureModel([Normal(mu,sigma), Exponential(lambda)], w)
    for i=1:length(d)
        d[i] ~ dist
    end
    n = length(d)
    return (d, w[1]*n, w[2]*n)
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