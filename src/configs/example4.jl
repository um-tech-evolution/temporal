# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example3
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example3
=#
export simtype
@everywhere simtype = 2    
@everywhere T = 2    # number of trials
@everywhere const N = 32        # Meta-population size
const num_subpops_list =[1,8]                     # Number of subpopulations
const num_attributes = 6        # number attributes for quantitative representation
const ngens = 12       # Generations after burn-in
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [1]                    # number emmigrants
#const mutation_stddev = 0.02
const mutation_stddev_list = [0.01]
const move_range=0.05
#const move_time_interval=10
const move_time_interval_list=[5]
#const horiz_select=false
const horiz_select_list=[true]
const uniform_start=true
const min_fit = 0.40
const linear_fitness=true
const linfit_slope_list=[1.0]
#const topology = "moore"
const topology_list = ["ring"]
const burn_in= 0.0    # generations of burn_in as a multiple of N
