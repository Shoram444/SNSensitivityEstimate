"""
    get_mean_conf(xs, chain, fit_function, alpha)

Used for plotting credible interval around the fit. 
The function calculates the hpdi for all samples per fit point. 
Returns three vectors:
   1. means: mean value of all samples at points xs
   2. lowers: mean - lower credible interval defined by alpha
   2. uppers: upper - mean credible interval defined by alpha
"""
function get_mean_conf(xs, chain, fit_function, alpha)
    
    n_pts = length(xs)
    n_samples= length(chain)

    # initiate an array that stores all the computed values
    # each row represents a sample of the posterior with n_pts as columns
    f_vals = Array{Float64}(undef, n_samples, n_pts )

    # initiate a parameters vector of size equal to number of parameters of the model
    # this will be filled at each iteration in for loop and used in fit_function
    n_params = length(chain.name_map.parameters)
    for s=1:n_samples
        params = [chain[s, i, 1] for i=1:n_params]
        f_vals[s,:] = [fit_function(x, params) for x in xs]
    end

    # for each column in the f_vals calculate the mean value, lower and upper interval
    m, l, u = Vector{Float64}(undef, n_pts), Vector{Float64}(undef, n_pts), Vector{Float64}(undef, n_pts)
    for (i,c) in enumerate(eachcol(f_vals))
        m_i = mean(c)
        l_i, u_i = MCMCChains._hpd(c; alpha= alpha)
        m[i] = m_i
        l[i] = m_i - l_i
        u[i] = u_i - m_i
    end
    return m, l, u
end


function get_median_conf(xs, chain, fit_function, alpha)
    
    n_pts = length(xs)
    n_samples= length(chain)

    # initiate an array that stores all the computed values
    # each row represents a sample of the posterior with n_pts as columns
    f_vals = Array{Float64}(undef, n_samples, n_pts )

    # initiate a parameters vector of size equal to number of parameters of the model
    # this will be filled at each iteration in for loop and used in fit_function
    n_params = length(chain.name_map.parameters)
    for s=1:n_samples
        params = [chain[s, i, 1] for i=1:n_params]
        f_vals[s,:] = [fit_function(x, params) for x in xs]
    end

    # for each column in the f_vals calculate the mean value, lower and upper interval
    m, l, u = Vector{Float64}(undef, n_pts), Vector{Float64}(undef, n_pts), Vector{Float64}(undef, n_pts)
    for (i,c) in enumerate(eachcol(f_vals))
        m_i = median(c)
        l_i, u_i = MCMCChains._hpd(c; alpha= alpha)
        m[i] = m_i
        l[i] = m_i - l_i
        u[i] = u_i - m_i
    end
    return m, l, u
end