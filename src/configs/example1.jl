# Configuration for running temporal simulation
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run.jl configs/example1
or
>  julia -p 4 -L TemporalEvolution.jl run.jl configs/example1
export simtype
@everywhere simtype = 2    
ii=#
@everywhere const N = 64        # Meta-population size
const num_subpops_list = [1,4]                     # Number of subpopulations
const num_attributes = 4        # number attributes for quantitative representation
const ngens = 1000       # Generations after burn-in
const ne = 1                    # number emmigrants
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [0,1]                    # number emmigrants
#const mutation_stddev = 0.02
const mutation_stddev_list = [0.01,0.02]
const move_range=0.1
#const move_time_interval=10
const move_time_interval_list=[10]
const opt_loss_cutoff=0.2
#const horiz_select=false
const horiz_select_list=[false,true]
const uniform_start=false
const burn_in= 1.0    # generations of burn_in as a multiple of N
