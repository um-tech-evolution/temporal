# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia run.jl configs/example3
or
>  julia -p 4 run.jl configs/example3
=#
simtype = 1    # run repeat_evolve_until_dead instead of repeat_evolve
num_trials = 4                    # number of trials with same parameter settings
N = 16        # Meta-population size
#num_subpops = [1,4]                     # Number of subpopulations
num_subpops = [1,4]                     # Number of subpopulations
num_attributes = 4        # number attributes for quantitative representation
#num_attributes = [4]        # number attributes for quantitative representation
ngens = 30       # Generations after burn-in
num_emigrants = 1                    # number emmigrants
#num_emigrants = [1]                    # number emmigrants
mutStddev = 0.01
#mutStddev = [0.01]
move_range=0.1
#move_time_interval=10
move_time_interval=[10]
probHSelect=1.0
horiz_select=false
#horiz_select=[false,true]
#horiz_select=[false]
#probHSelect=1.0
probHSelect=[1.0]
#horiz_mutate=true
horiz_mutate=[true]
uniform_start=true
minFit = 0.40
ideal_init=0.5
#linear_fitness=true
linfit_slope = 1.0
topology= ["ring"]
burn_in= 0.0    # generations of burn_in as a multiple of N
additive_error=true
