# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example1
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example1
export simtype
=#
simtype = 2    # run  repeat_evolve instead of repeat_evolve_until_dead
num_trials = 1  # number of trials
N = 8        # Meta-population size
num_subpops = [1,4]                     # Number of subpopulations
#num_attributes = 4        # number attributes for quantitative representation
num_attributes = [4]        # number attributes for quantitative representation
ngens = 2       # Generations after burn-in
#num_emigrants = 1                    # number emmigrants
num_emigrants = [1]                    # number emmigrants
#mutStddev = 0.03
mutStddev = [0.03]
move_range=0.1
move_time_interval=10
#move_time_interval=[10]
horiz_select=false
#horiz_select=[false]
#probHSelect=1.0
probHSelect=[1.0]
horiz_mutate=[false]
uniform_start=true
ideal_init=0.5
minFit = 0.25
linfit_slope=1.0
#linfit_slope=[1.0]
topology= ["ring"]
burn_in= 0.2    # generations of burn_in as a multiple of N
additive_error=true
