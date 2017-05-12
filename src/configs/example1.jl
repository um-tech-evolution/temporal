# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example1
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example1
export simtype
=#
@everywhere simtype = 2    
const T = 1
@everywhere const N = 27        # Meta-population size
const num_subpops_list = [9]                     # Number of subpopulations
#const num_attributes = 2        # number attributes for quantitative representation
const num_attributes_list = [2]        # number attributes for quantitative representation
const ngens = 6       # Generations after burn-in
const ne = 1                    # number emmigrants
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [1]                    # number emmigrants
#const mutStddev = 0.02
const mutStddev_list = [0.01]
const move_range=0.1
#const move_time_interval=10
const move_time_interval_list=[10]
const opt_loss_cutoff=0.2
const probHSelect=0.2
#const horiz_select=false
const horiz_select_list=[true]
const uniform_start=false
const minFit = 0.0
const topology_list=["moore"]
const linear_fitness=false
const linfit_slope=1.0
#const linfit_slope_list=[1.0]
const burn_in= 1.0    # generations of burn_in as a multiple of N
