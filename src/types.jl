export variant_type, fitness_location_type, Population, PopList, temporal_result_type, subpop_properties_type, variant_type, fitness_location_type,
    subpop_properties
using Distributions
using StatsBase
typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

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
  mutation_stddev::Float64  # standard deviation of normal distribution of mutation perturbations
  ideal_init::Float64     # initial values of ideal corresponding to initial optima
  move_range::Float64     # on optima change, a uniform random number in the interval [-move_range,move_range] is added to each ideal
  move_time_interval::Int64  # The optima are moved every move_time_interval generations
  min_fit::Float64      # Lower bound on fitness
  linear_fitness::Bool      # if true, fitness is 0.5 - Euclidean distance
  linfit_slope::Float64
  topology::String      # Neighborhood topology for horizontal_transfer_by_fitness: must be "circular", "ring", "vonneumann", "moore", "random", or "none"
  mean_fraction_subpops_below_min_fit::Float64  # mean over generations of the fraction of subpops with an individual whose fitnsss is greater than min_fit
  fraction_gens_with_all_subpops_below_min_fit::Float64
  fitness_mean::Float64
  fitness_variance::Float64
  attribute_variance::Float64
  generational_lifetime::Float64
  move_update_lifetime::Float64
  gen_limit_reached_count::Int64
end

