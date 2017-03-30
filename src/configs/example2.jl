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
const num_subpops_list = [1,4,8]                     # Number of subpopulations
const num_attributes = 4        # number attributes for quantitative representation
const ngens = 100       # Generations after burn-in
#const num_emmigrants = 1                    # number emmigrants
const num_emmigrants_list = [0,1,2]                    # number emmigrants
const normal_stddev = 0.01
const move_range=0.1
const move_time_interval=10
const opt_loss_cutoff=0.2
#const horiz_select=false
const horiz_select_list=[false]
const uniform_start=false
const burn_in= 0.1    # generations of burn_in as a multiple of N
