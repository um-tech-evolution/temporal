# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia run.jl examples/ed_example1
or
>  julia -p 4 run.jl examples/ed_example1
=#
simtype = 1    # run repeat_evolve_until_dead instead of repeat_evolve
num_trials = 2                    # number of trials with same parameter settings
N = 32        # Meta-population size
num_subpops = 8                     # Number of subpopulations
num_attributes = 4        # number attributes for quantitative representation
ngens = 5       # Generations after burn-in
num_emigrants = 1                    # number emmigrants
migration_rate = 0.0
mutStddev = 0.05
move_range=0.1
move_time_interval=3
probHSelect=1.0
horiz_select=true
probHSelect=1.0
horiz_mutate=true
uniform_start=true
minFit = 0.40
ideal_init = 0.5
#linear_fitness=true
linfit_slope = 1.0
topology = "ring"
burn_in= 0.0    # generations of burn_in as a multiple of N
additive_error=true
