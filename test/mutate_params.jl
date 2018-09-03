# Used by test_horiz.jl
# Configuration for to run unit tests for functions in horiz.jl
#=
Recommended command line to run:
>  julia test_horiz.jl 
=#
simtype = 1    # run repeat_evolve_until_dead instead of repeat_evolve
num_trials = 1                    # number of trials with same parameter settings
probHSelect=1.0
N = 8         # Meta-population size
num_subpops = 4                     # Number of subpopulations
num_attributes = 2        # number attributes for quantitative representation
ngens = 20       # Generations after burn-in
num_emigrants = 1                    # number emmigrants
migration_rate = 0.0
mutStddev = 0.10
move_range=0.04
move_time_interval=2
horiz_select=false
horiz_mutate=false
uniform_start=false
minFit = 0.40
ideal_init = 0.5
#linear_fitness=true
linfit_slope = 1.0
topology = "global"
burn_in= 0.0    # generations of burn_in as a multiple of N
additive_error=true
