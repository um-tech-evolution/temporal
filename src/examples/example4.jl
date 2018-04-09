# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example3
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example3
=#
export simtype
@everywhere simtype = 2    # run  repeat_evolve instead of repeat_evolve_until_dead
@everywhere T = 2    # number of trials
@everywhere const N = 16       # Meta-population size
const num_subpops_list =[4]                     # Number of subpopulations
#const num_attributes = 6        # number attributes for quantitative representation
const num_attributes_list = [3]        # number attributes for quantitative representation
const ngens = 4       # Generations after burn-in
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [1]                    # number emmigrants
#const mutStddev = 0.02
const mutStddev_list = [0.01]
const move_range=0.05
#const move_time_interval=10
const move_time_interval_list=[5]
const horiz_select=false
#const horiz_select_list=[false]
#const probHSelect=1.0
const probHSelect_list=[1.0]
const uniform_start=true
const minFit = 0.40
const linear_fitness=true
#const linfit_slope_list=[1.0]
const linfit_slope=1.0
#const topology = "moore"
const topology_list = ["global"]
const burn_in= 0.0    # generations of burn_in as a multiple of N
