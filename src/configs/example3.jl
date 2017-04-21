# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example3
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example3
=#
export simtype
@everywhere simtype = 2    
const T = 1                    # number of trials with same parameter settings
@everywhere const N = 64        # Meta-population size
#const num_subpops_list = [1,8,16]                     # Number of subpopulations
const num_subpops_list = [1,16]                     # Number of subpopulations
const num_attributes = 4        # number attributes for quantitative representation
const ngens = 10       # Generations after burn-in
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [0,1]                    # number emmigrants
#const mutation_stddev = 0.02
const mutation_stddev_list = [0.01]
const move_range=0.1
#const move_time_interval=10
const move_time_interval_list=[10]
#const horiz_select=false
#const horiz_select_list=[false,true]
const horiz_select_list=[false]
const uniform_start=true
const min_fit = 0.35
const linear_fitness=true
const linfit_slope_list = [0.75]
const topology_list= ["none","circular","ring","vonneumann","moore","global"]
const burn_in= 0.0    # generations of burn_in as a multiple of N
