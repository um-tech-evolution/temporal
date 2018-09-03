# Configuration for running spatial simulation
#=
Recommended command line to run:
>  julia -L run.jl configs/example1
or
>  julia -p 4 run.jl configs/example1
=#
simtype = 4    
num_trials = 2
N = 8        # Meta-population size
num_subpops = [2]                     # Number of subpopulations
mu = 0.00                 # per-individual innovation rate 
num_emigrants = 1                    # number emigrants
num_attributes = 1        # number attributes for quantitative representation
use_fit_locations=false
horiz_select=false
extreme_variation=false
circular_variation=true
ngens = 2       # Generations after burn-in
burn_in= 3    # if integer, generations of burn in.  If float, generations of burn_in as a multiple of N
mutStddev = 0.04
ideal_max = 0.8
ideal_min = 0.2
ideal_range = 0.1
linfit_slope = 1.0
additive_error = false
neutral = false
#num_fit_locations=0
