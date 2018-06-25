export variant_type, fitness_location_type, Population, PopList, temporal_result_type, subpop_properties_type, variant_type, 
    fitness_location_type, subpop_properties, generational_innovation_counts
using Distributions
using StatsBase
const Population =  Array{Int64,1}
const PopList = Array{Population,1}

const param_type = Dict{Symbol,Any}
const result_type = Dict{Symbol,Any}

type variant_type
  fitness::Float64    # The fitness of this variant
  attributes::Vector{Float64}   # attributes of the variant
end

type fitness_location_type
  ideal::Vector{Float64}   # Ideal values for attributes in the environment of this subpop
end

type subpop_properties
  generational_subpop_alive::Vector{Bool}
  prev_subpop_alive::Vector{Bool}
  current_subpop_alive::Vector{Bool}
end

type generational_innovation_counts
  pos::Int64
  half_pos::Int64
  neg::Int64
  half_neg::Int64
end

type trial_innovation_counts
  pos::Float64
  half_pos::Float64
  neg::Float64
  half_neg::Float64
end

global temporal_param_fields = [
  :simtype    # if simtype==2, run  repeat_evolve. if symtype==1, run repeat_evolve_until_dead for length of retention.
  :simname     # name of configuration file (with .jl extension) and of output file (with .csv extenseion)
  :num_trials  # number of trials
  :N        # Meta-population size
  :num_subpops                     # Number of subpopulations
  :num_attributes         # number attributes for quantitative representation
  :ngens                  # Generations after burn-in
  :num_emigrants          # number emigrants in horizontal transmission between subpopulations, if nonzero, migration_rate must be 0
  :migration_rate         # Float fraction of subpop to emigrate.  If nonzero, num_emigrants must be 0
  :mutStddev              # Float.  Standard deviation of normal random deviate of copy error
  :move_range             # Float.  Add a uniform random number to each attribute ideal to move optima
  :move_time_interval  # Integer.  Number of generations between moving ideals
  :horiz_select      # Boolean.  If true, use proportional selection to select members of source pop to emigrate
  :horiz_mutate      # Boolean.  If true, mutate attributes as part of horizontal transmission
  :probHSelect       # Float.  Probability of selecting the best subpopulation as the source of emigration
  :uniform_start
  :ideal_init  # should be 0.5 for additive error, 1.0 for multiplicative error
  :minFit      # Float.  minimum fitness.  When an individual has this fitness, it has "lost" the ideal
  :linfit_slope  # 1.0 for most runs
  :topology   # one of circular, ring, moore, vonneumann, global.  Circular is not by fitness, the others are.
  :burn_in    # generations of burn_in as a multiple of N
  :additive_error   # if additive_error==false, multiplicative error
]

global temporal_result_fields = [
  :fitness_mean
  :fitness_variance
  :attribute_variance
  :stddev_fitness
  :stddev_attributes
  :mean_fraction_subpops_below_minFit
  :fraction_gens_with_all_subpops_below_minFit
  :innovations_per_gen_trial
  :half_innovations_per_gen_trial
  :deleterious_mutations_per_gen_trial
  :half_deleterious_mutations_per_gen_trial
]

global subpop_alive_result_fields = [
  :generational_lifetime
  :move_update_lifetime
  :gen_limit_reached_count
]


const temporal_result_type = Dict{Symbol,Any}
