# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example1
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example1
=#
simtype = 2   # run  repeat_evolve instead of repeat_evolve_until_dead
num_trials = 3
N = 4        # Meta-population size
num_subpops=[1]                     # Number of subpopulations
num_attributes = 1        # number attributes for quantitative representation
ngens = 3       # Generations after burn-in
num_emigrants = 0                    # number emigrants
mutStddev = 0.01
move_range=0.1
move_time_interval=0
probHSelect=[0.0]
horiz_select=true
horiz_mutate=true
uniform_start=false
ideal_init=0.5
minFit = 0.0
topology=["vonneumann"]
linfit_slope=[1.0]
burn_in= 0.0    # generations of burn_in as a multiple of N
