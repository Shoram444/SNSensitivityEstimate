const ALL_PROCESSES = vcat(signal, background)
const NPARAM = 2 * length(VARS_)   # center + logwidth per var

const ROI_BUFFERS = [zeros(Float64, NPARAM_) for _ in 1:Threads.nthreads()]

function normalized_to_roi_vector!(roi, x̂)

    WMIN = 0.05   # minimum fractional width

    @inbounds for (i, var) in enumerate(VARS_)

        L, U = BOUNDS_[Symbol(var)]
        R = U - L

        center_frac = clamp(x̂[2i-1], 0.0, 1.0)
        logw_frac   = clamp(x̂[2i],   0.0, 1.0)

        # log-uniform width sampling
        wmin = WMIN
        wmax = 1.0

        logw = log(wmin) + logw_frac*(log(wmax) - log(wmin))
        w = exp(logw)

        center = L + center_frac*R
        halfw = 0.5*w*R

        minv = clamp(center - halfw, L, U)
        maxv = clamp(center + halfw, L, U)

        # enforce strict ordering
        if maxv <= minv
            maxv = min(minv + 1e-12, U)
        end

        roi[2i-1] = minv
        roi[2i]   = maxv
    end
end

function prob_fast(x::AbstractVector)

    roi = ROI_BUFFERS[Threads.threadid()]

    normalized_to_roi_vector!(roi, x)

    val = -get_s_to_b(
        SNparams,
        α,
        ALL_PROCESSES,
        roi;
        approximate="formula"
    )

    return isfinite(val) ? val : 1e20
end

function prob_parallel(X::AbstractMatrix)

    fitness = Vector{Float64}(undef, size(X,1))

    Threads.@threads for i in axes(X,1)
        fitness[i] = prob_fast(view(X,i,:))
    end

    return fitness
end

LB = zeros(NPARAM)
UB = ones(NPARAM)

bounds = boxconstraints(lb=LB, ub=UB)

LB = zeros(NPARAM)
UB = ones(NPARAM)

bounds = boxconstraints(lb=LB, ub=UB)

options = Options(
    parallel_evaluation = true,
    verbose = true,
    time_limit = 120.0
)

algo = ECA(N=30, K=5; options)

x0 = reshape(repeat([0.5, 1.0], length(VARS_)), 1, :)

set_user_solutions!(algo, x0, prob_fast)

result = optimize(prob_parallel, bounds, algo)
