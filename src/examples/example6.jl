# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia run.jl configs/example3
or
>  julia -p 4 run.jl configs/example3
=#
simtype = 1    # run repeat_evolve_until_dead instead of repeat_evolve
num_trials = 2                    # number of trials with same parameter settings
N = 8        # Meta-population size
num_subpops = [4]                     # Number of subpopulations
num_attributes = 2        # number attributes for quantitative representation
#num_attributes_list = [2]        # number attributes for quantitative representation
ngens = 10       # Generations after burn-in
num_emigrants = 1                    # number emmigrants
#num_emigrants_list = [1]                    # number emmigrants
mutStddev = 0.01
#mutStddev_list = [0.01]
move_range=0.1
move_time_interval=3
#move_time_interval_list=[3]
horiz_select=false
#horiz_select_list=[false]
probHSelect=0.0
#probHSelect_list=[0.0]
horiz_mutate=false
#horiz_mutate_list=[false]
uniform_start=true
ideal_init = 0.5
minFit = 0.40
#linear_fitness=true
linfit_slope = 1.0
topology= ["ring"]
burn_in= 0.5    # generations of burn_in as a multiple of N
additive_error=true
