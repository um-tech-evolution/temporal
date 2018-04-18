# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example1
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example1
export simtype
=#
@everywhere simtype = 2    # run  repeat_evolve instead of repeat_evolve_until_dead
const T = 1  # number of trials
@everywhere const N = 8        # Meta-population size
const num_subpops_list = [1,4]                     # Number of subpopulations
#const num_attributes = 4        # number attributes for quantitative representation
const num_attributes_list = [4]        # number attributes for quantitative representation
const ngens = 2       # Generations after burn-in
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [1]                    # number emmigrants
#const mutStddev = 0.01
const mutStddev_list = [0.03]
const move_range=0.1
#const move_time_interval=10
const move_time_interval_list=[10]
const horiz_select=false
#const horiz_select_list=[false]
#const probHSelect=1.0
const probHSelect_list=[1.0]
const horiz_mutate_list=[false]
const uniform_start=true
const minFit = 0.25
const linfit_slope=1.0
#const linfit_slope_list=[1.0]
const topology_list= ["ring"]
const burn_in= 0.2    # generations of burn_in as a multiple of N

