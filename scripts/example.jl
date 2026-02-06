#############################################
#############################################
##
##  LOAD DEPENDENCIES
##
#############################################
#############################################
using SNSensitivityEstimate, UnROOT, Metaheuristics

const α= 1.64485362695147 # 90% C.L. constant

#############################################
## LOAD SIMULATION Data 
#############################################

# here of course put your own path, maybe absolute rather than relative
signal_path = "data/sims/final_phd/fal5_12perc_Boff_Cimrman_J41/bb0nu_foil_bulk.root"
tree_name = "tree"
f1 = ROOTFile(signal_path) 
d1 = LazyTree(f1, tree_name, keys(f1[tree_name])) 

bkg1_path = "data/sims/final_phd/fal5_12perc_Boff_Cimrman_J41/bb_foil_bulk.root"
f2 = ROOTFile(bkg1_path)
d2 = LazyTree(f2, tree_name, keys(f2[tree_name]))

bkg2_path = "data/sims/final_phd/fal5_12perc_Boff_Cimrman_J41/Bi214_wire_surface.root"
f3 = ROOTFile(bkg2_path)
d3 = LazyTree(f3, tree_name, keys(f3[tree_name]))

# and so on, add all backgrounds you want to consider

#############################################
## INITIATE VARIABLES TO READ AND RANGES
#############################################
vars = [
    "phi", 
    "sumE", 
    # Here you should add at least dy, dz, deltaCalo
]

var_range = (
    phi = (0,180),
    sumE = (0, 3500),
)

##############################################
## INITIATE PROCESSES
##############################################


signal_process = DataProcessND(
    d1, # data
    "nu0bb", # isotopeName
    true, # signal
    1.0, # acitivity 
    1.0 * 365 * 24 * 3600, # time in seconds
    1e7, # nsimulated events
    bins, # bins 
    6.11, # amount
    vars # variable names
)

background_process1 = DataProcessND(
    d2,
    "bb2nu",
    false,
    100/1e3,
    1.0 * 365 * 24 * 3600,
    1e8,
    bins,
    6.11,
    vars
)

background_process2 = DataProcessND(
    d4,
    "radon",
    false,
    150/2e6, # 150 uBq/m3
    1.0 * 365 * 24 * 3600, # 1 year of measurement
    1e8, # 1e8 simulated events 
    bins,
    15, # 15 m3 of tracking volume
    vars
)


##################################
#### PROCESSESS
##################################

# list all processes here and combine them in common vector
all_processes = [
    signal_process,
    background_process1,
    background_process2,
]

# Test simple roi, you can play around by hand to get the feeling
my_roi = (
    phi = (0, 140),
    sumE = (2000, 3100),
)

# there are 2 approximation methods "table" and "formula", 
# table is more precise but takes longer to calculate, 
# formula is fine for quick tests
my_sensitivity = get_sensitivityND(SNparams, α, all_processes, my_roi; approximate="table") 

# print your result, it shows the sensitivity, efficiency, background and roi
print(my_sensitivity) 



##############################################
## OPTIMIZATION BY METAHEURISTICS.jl
##############################################

# First we have to create a function that takes a proposed roi and returns the signal over background ratio -- we call this the "prob" function
# this is the function that we want to maximize, so we return the negative value (because Metaheuristics.jl minimizes functions)
# the proposed_roi is a vector of values that we have to convert to a roi tuple
prob(proposed_roi) = - get_s_to_b(SNparams, α, all_processes, proposed_roi; approximate="formula")


# Next we have to define the search ranges for each variable, 
# this is done automatically by providing a process (can be any process, here we use the signal process)
# the output is a 2*n vector of ranges, where n is the number of variables -- we have 2 values per variable: lower and upper bound
searchRange = make_stepRange(signal_process)


# Next we define from the searchRange the lower and upper bounds for the optimization
lower_bound = [x[1] for x in searchRange] .|> float
upper_bound = [x[2] for x in searchRange] .|> float


##############################################
# Now we can set up the optimization problem and run it
##############################################

# first the optimization options, basically the hyperparameters of the optimization algorithm
# you can leave them empty - default values will be used,
# or you can try to tune yourself to see what the results give, basically you want to find a good balance between time and quality of results
# For sure set up at least a time_limit and verbosity to true to see progress
# The full list of parameters can be found in https://jmejia8.github.io/Metaheuristics.jl/stable/api/#Metaheuristics.Options
# Use chatGPT to help you if needed
options = Options(;
    # x_tol = 1.0, # tolerance in position space
    # f_tol = 1e-5, # tolerance in function value
    # f_tol_rel = 1e-5, # relative tolerance in function value
    # f_tol_abs = 1e-5, # absolute tolerance in function value
    time_limit = 60*3.0, # time limit in seconds, must include .0 to be a Float64 value, doesn't take integers
    # parallel_evaluation = true, # parallel evaluation is cool, but needs more setup and only works with few optimizers
    verbose = true, # verbosity on
    # iterations = 15, # number of iterations / generations
    # store_convergence = true # store convergence history
)

# Set up boxconstraints bounds - just the lower and upper bounds we defined before
bounds = boxconstraints(lb = lower_bound, ub = upper_bound)

# initial guess for the roi, this can be left out, but can help the optimizer to start from a reasonable point
x0 = float.([
    5, 180, # phi
    2700, 3000, # sumE
    # of course add all the variables you defined earlier
])


# Set up optimizer algorithm, I generally use ECA - supports parallel evaluation, PSO - also can use parallel, SA - single core, but pretty good
# The algorithm is defined by calling the function with the options as argument, if you want ECA you'd do algo = ECA(;options), etc.
algo = PSO(;options)

# Set the initial guess for the optimizer
set_user_solutions!(algo, x0, prob)


# Now we run the optimization, depending on how much time you gave it in the options, this can take a while
result = optimize(prob, bounds, algo)

# We can show the best result found
@show minimum(result)

# and store it's value as the best roi
@show res=  minimizer(result)


# then we conver the result to roi tuple and calculate the sensitivity for it
best = get_best_ROI_ND(res, signal)
best_sens = get_sensitivityND(SNparams, α, vcat(signal, background), best; approximate="table")

# finally print the result
print(best_sens)