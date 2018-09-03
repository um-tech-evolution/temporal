# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia run.jl examples/example2
or
>  julia -p 4 run.jl examples/example2
=#
# simtype=4 means spatial structure by changing the ideal values for attributes
simtype = 4
num_trials = 4
N = 8        # Meta-population size
#num_subpops = [1,2]
num_subpops = 1
mutStddev = 0.1
#mu = 0.05                 # per-individual innovation rate 
mu = 0.0                 # per-individual innovation rate 
#num_emigrants = [0,2]                    # number emigrants
num_emigrants = 0                    # number emigrants
num_attributes = 2        # number attributes for quantitative representation
ngens = 2       # Generations after burn-in
#horiz_select=[false,true]
horiz_select=false
circular_variation=true
extreme_variation=false
burn_in= 0    # generations of burn_in as a multiple of N
#use_fit_locations=[false,true]
use_fit_locations=false
ideal_max = 0.8
ideal_min = 0.2
ideal_range = 0.1
linfit_slope = 0.0
neutral = false
additive_error=true
