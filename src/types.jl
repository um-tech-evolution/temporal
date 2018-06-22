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
  :simtype    # run  repeat_evolve instead of repeat_evolve_until_dead
  :simname     # name of configuration file (with .jl extension) and of output file (with .csv extenseion)
  :num_trials  # number of trials
  :N        # Meta-population size
  :num_subpops                     # Number of subpopulations
  :num_attributes         # number attributes for quantitative representation
  :ngens                  # Generations after burn-in
  :num_emigrants          # number emigrants in horizontal transmission between subpopulations
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
  :topology   # one of circular, ring, moore, vonneumann, global
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
#=
type temporal_result_type
  simtype::Int64   # 1:  evolve_until_dead,   2: evolve with fixed generations
  num_trials::Int64   # number of repeated runs with the same parameter values
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  ne::Int64  # number of emmigrants in horizontal transfer
  num_attributes::Int64  # number of attributes of a variant
  ngens::Int64  # number of generations after burn-in
  burn_in::Float64
  uniform_start::Bool      # if true, initial population consists of copies of ideal, otherwise random
  horiz_select::Bool       # Whether to use selection during horzontal transfer
  probHSelect::Float64     # Whether to choose higher fitness source for horiz trans or random source
  mutStddev::Float64  # standard deviation of normal distribution of mutation perturbations
  ideal_init::Float64     # initial values of ideal corresponding to initial optima
  move_range::Float64     # on optima change, a uniform random number in the interval [-move_range,move_range] is added to each ideal
  move_time_interval::Int64  # The optima are moved every move_time_interval generations
  minFit::Float64      # Lower bound on fitness
  linfit_slope::Float64    # the slope of the peak for linear fitness.  Larger means a smaller, steeper, peak
  horiz_mutate::Bool      # if true, do mutation (copy error) is done as part of horizontal transmission
  topology::String      # Neighborhood topology for horizontal_transfer_by_fitness: must be "circular", "ring", "vonneumann", "moore", "random", or "none"
  mean_fraction_subpops_below_minFit::Float64  # mean over generations of the fraction of subpops with an individual whose fitnsss is greater than minFit
  fraction_gens_with_all_subpops_below_minFit::Float64
  fitness_mean::Float64
  fitness_variance::Float64
  attribute_variance::Float64
  generational_lifetime::Float64
  move_update_lifetime::Float64
  gen_limit_reached_count::Int64
  innovations_per_gen_trial::Float64
  half_innovations_per_gen_trial::Float64
  deleterious_mutations_per_gen_trial::Float64
  half_deleterious_mutations_per_gen_trial::Float64
end

function temporal_result( simtype::Int64, num_trials::Int64, N::Int64, num_attributes::Int64, num_subpops::Int64, ngens::Int64, mutStddev::Float64, 
    num_emmigrants::Int64, move_range::Float64, move_time_interval::Int64, horiz_select::Bool=false, probHSelect::Float64=1.0, minFit::Float64=0.0; 
    horiz_mutate::Bool=false, topology::String="circular", 
    uniform_start::Bool=false, linfit_slope::Float64=1.0, burn_in::Float64=1.0 )
  ideal_init = 0.5
  return temporal_result_type( simtype, num_trials, N, num_subpops, num_emmigrants, num_attributes, ngens, burn_in, uniform_start, horiz_select, probHSelect,
      mutStddev, ideal_init, move_range, move_time_interval, minFit,  linfit_slope, horiz_mutate, topology, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0 )
end
=#
