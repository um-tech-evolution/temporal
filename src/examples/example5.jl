# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example3
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example3
=#
simtype = 1    # run repeat_evolve_until_dead instead of repeat_evolve
T = 4                    # number of trials with same parameter settings
num_trials = 4                    # number of trials with same parameter settings
N = 64        # Meta-population size
num_subpops_list = [16]                     # Number of subpopulations
num_subpops = [16]                     # Number of subpopulations
num_attributes = [4]        # number attributes for quantitative representation
num_attributes_list = [4]        # number attributes for quantitative representation
ngens = 30       # Generations after burn-in
num_emmigrants = [1]                    # number emmigrants
num_emmigrants_list = [1]                    # number emmigrants
mutStddev = 0.01
mutStddev_list = [0.01]
move_range=0.1
move_time_interval=[10]
move_time_interval_list=[10]
probHSelect=1.0
horiz_select=true
#horiz_select_list=[true]
probHSelect=1.0
probHSelect_list=[1.0]
horiz_mutate=true
horiz_mutate_list=[true]
uniform_start=true
minFit = 0.40
ideal_init = 0.5
linear_fitness=true
#linfit_slope_list = [0.8]
linfit_slope = 1.0
topology = ["ring"]
topology_list= ["ring"]
burn_in= 0.0    # generations of burn_in as a multiple of N
