# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia run.jl configs/example3
or
>  julia -p 4 run.jl configs/example3
=#
simtype = 2    # run  repeat_evolve instead of repeat_evolve_until_dead
num_trials = 1                    # number of trials with same parameter settings
N = 64        # Meta-population size
#num_subpops = [1,8,16]                     # Number of subpopulations
num_subpops = [1,16]                     # Number of subpopulations
#num_attributes = 4        # number attributes for quantitative representation
num_attributes = [4]        # number attributes for quantitative representation
ngens = 10       # Generations after burn-in
#num_emmigrants = 1                    # number emmigrants
num_emigrants = [0,1]                    # number emmigrants
#mutStddev = 0.02
mutStddev = [0.01]
move_range=0.1
#move_time_interval=10
move_time_interval=[10]
horiz_mutate=false
horiz_select=false
#horiz_select=[false,true]
#probHSelect=1.0
probHSelect=[0.5,1.0]
uniform_start=true
minFit = 0.35
ideal_init = 0.5
#linear_fitness=true
#linfit_slope = [0.75]
linfit_slope = 1.0
topology= ["none","circular","ring","vonneumann","moore","global"]
burn_in= 0.0    # generations of burn_in as a multiple of N
additive_error=true
