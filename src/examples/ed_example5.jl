# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia run.jl examples/ed_example5
or
>  julia -p 4 run.jl examples/ed_example5
=#
simtype = 1    # run repeat_evolve_until_dead instead of repeat_evolve
num_trials = 6                    # number of trials with same parameter settings
N = [36,64,100]        # Meta-population size
num_subpops = "sqrt"                     # Number of subpopulations
num_attributes = 1        # number attributes for quantitative representation
ngens = 1000       # Generations after burn-in
num_emigrants = 0                    # number emmigrants
migration_rate = 0.1
mutStddev = 0.10
move_range=0.04
move_time_interval=2
probHSelect=0.0
horiz_select=true
probHSelect=0.0
horiz_mutate=true
uniform_start=true
minFit = 0.40
ideal_init = 0.5
#linear_fitness=true
linfit_slope = 1.0
topology = "ring"
burn_in= 0.0    # generations of burn_in as a multiple of N
additive_error=false
