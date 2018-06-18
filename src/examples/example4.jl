# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia run.jl configs/example3
or
>  julia -p 4 run.jl configs/example3
=#
simtype = 2    # run  repeat_evolve instead of repeat_evolve_until_dead
num_trials = 2    # number of trials
N = 16       # Meta-population size
num_subpops =[4]                     # Number of subpopulations
#num_attributes = 6        # number attributes for quantitative representation
num_attributes = [3]        # number attributes for quantitative representation
ngens = 4       # Generations after burn-in
#num_emigrants = 1                    # number emmigrants
num_emigrants = [1]                    # number emmigrants
#mutStddev = 0.02
mutStddev = [0.01]
move_range=0.05
#move_time_interval=10
move_time_interval=[5]
horiz_mutate=[false,true]
horiz_select=false
#horiz_select=[false]
#probHSelect=1.0
probHSelect=[1.0]
uniform_start=true
minFit = 0.40
ideal_init = 0.5
#linear_fitness=true
#linfit_slope=[1.0]
linfit_slope=1.0
#topology = "moore"
topology = ["global"]
burn_in= 0.0    # generations of burn_in as a multiple of N
additive_error=true
