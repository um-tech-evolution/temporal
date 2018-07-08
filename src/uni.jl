# Legal fields for the dictionaries that are used to parameters, results, and individuals.
# Suggested test sequence:
# julia> paramd = init_dictionary( param_fields )   # builds dictionary paramd
# julia> read_parameter_file("param",paramd)        # adds parameters from param file into dictionary
# julia> build_pmap_list( paramd )                  # builds list of dictionaries for each trial

export temporal_param_fields, init_dictionary, read_parameter_file, print_meta_pop, build_pmap_list, print_meta_pop_attributes

# Parameter fields for temporal 
#=  moved to types.jl
const temporal_param_fields = [ 
  :simtype,   # 1:  evolve_until_dead,   2: evolve with fixed generations 
  :num_trials,   # number of repeated runs with the same parameter values 
  :N,   # meta-population size 
  :num_subpops,   # number of subpopulations 
  :ne,  # number of emigrants in horizontal transfer 
  :num_emigrants,  # number of emigrants
  :num_attributes,  # number of attributes of a variant 
  :ngens,  # number of generations after burn-in
  :burn_in,
  :uniform_start,      # if true, initial population consists of copies of ideal, otherwise random
  :horiz_select,       # Whether to use selection during horzontal transfer
  :probHSelect,     # Whether to choose higher fitness source for horiz trans or random source
  :mutStddev,  # standard deviation of normal distribution of mutation perturbations
  :ideal_init,     # initial values of ideal corresponding to initial optima
  :move_range,     # on optima change, a uniform random number in the interval [-move_range,move_range] is added to each ideal
  :move_time_interval,  # The optima are moved every move_time_interval generations
  :minFit,      # Lower bound on fitness
  :linfit_slope,    # the slope of the peak for linear fitness.  Larger means a smaller, steeper, peak
  :horiz_mutate,      # if true, do mutation (copy error) is done as part of horizontal transmission
  :topology      # Neighborhood topology for horizontal_transfer_by_fitness: must be "circular", "ring", "vonneumann", "moore", "random", or "none"
]
=#

# param fields for testing
const param_fields = [
  :N,
  :minFit,
  :num_trials
]

# Results for temporal
const result_fields = [
  :mean_fraction_subpops_below_minFit,  # mean over generations of the fraction of subpops with an individual whose fitnsss is greater than minFit
  :fraction_gens_with_all_subpops_below_minFit,
  :fitness_mean,
  :fitness_variance,
  :attribute_variance,
  :generational_lifetime,
  :move_update_lifetime
]

# Returns a parameter or result dictionary with the specified fields set with value :null
# Example:  paramd = init_dictionary( param_fields )
# Example:  resultd = init_dictionary( result_fields )
function init_dictionary( field_list )
  dict = Dict{Symbol,Any}()
  for f in field_list
    dict[f] = :null
  end
  dict
end

function print_meta_pop( paramd::param_type, meta_pop::PopList, vt::Dict{Int64,variant_type} )
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  for  j = 1:paramd[:num_subpops]
    print(" ")
    #=
    for k in meta_pop[j]
      @printf("%2d %5.3f ",k,vt[k].fitness)
    end
    =#
    for k in meta_pop[j]
      above_minFit = vt[k].fitness > paramd[:minFit] ? 1 : 0
      print( above_minFit )
    end
    #print(subpop_alive(meta_pop[j],vt,paramd[:minFit]),"] ")
  end
  println()
end

function print_meta_pop_attributes( paramd::param_type, meta_pop::PopList, vt::Dict{Int64,variant_type} )
  println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  for  j = 1:paramd[:num_subpops]
    print("[")
    for k in meta_pop[j]
      @printf("[%d %5.3f ",k,vt[k].fitness)
      attr = vt[k].attributes
      print(attr)
      #=
      for a in attr
        @printf("%5.3f ",a)
      end
      =#
      print("]")
    end
    println("]")
    #print(subpop_alive(meta_pop[j],vt,paramd[:minFit]),"] ")
  end
  println()
end


# Returns the dictionary param_dict with the parameter values specified in the parameter file inserted into the param_dictionary.
# The param_dictionary param_dict must have been initialized with all of the parameters specified in the parameter file.
function read_parameter_file( filename::AbstractString, param_dict::Dict{Symbol,Any} )
  open(filename,"r") do f
    lines = readlines(f)
    comment_line = false
    for line in lines
      #println("line: ",line)
      if comment_line && strip(line)[1:2] != "=#"
        continue
      elseif comment_line
        comment_line = false
        continue
      end
      if strip(line)[1:2] == "#="
        comment_line = true
      end
      if line[1] != '#' 
        pos_equals = searchindex(line,"=")   # position of the = symbol
        #println("pos_equals: ",pos_equals)
        pos_comment_char = searchindex(line,"#")  # position of a comment char after the value assignment
        #println("pos_comment_char: ",pos_comment_char)
        if pos_comment_char == 0
          value = parse(strip(line[pos_equals+1:end]))
          #println("n value: ",value)
        else
          value = parse(strip(line[pos_equals+1:pos_comment_char-1]))
          #println("c value: ",value)
        end
        #println("value: ",value)
        if typeof(value) == Expr
          value = eval(value)
        end
        key = parse(strip(line[1:pos_equals-1]))
        #println("key: ",key,"  value: ",value)
        if !haskey(param_dict,key)
          error("dictionary not initialized for key: ",key)
        end
        param_dict[key] = value
      end
    end
  end
  for key in keys(param_dict)
    if param_dict[key] == :null
      println("warning:  parameter: ",key," not initialized")
    end
  end
  param_dict
end  

# builds a list of dictionaries that can be used to run a sequence of trials with different
#   parameter settings using pmap()  (parallel map).
function build_pmap_list(  param_dict::Dict{Symbol,Any} )
  array_field_list = Symbol[]
  for k in keys(param_dict)
    if typeof(param_dict[k]) <: Array
      Base.push!( array_field_list, k )
    end
  end
  #println(array_field_list)
  result = []
  n = length(array_field_list)
  plist_lengths = [ length( param_dict[k] ) for k in array_field_list ]
  #println("plist_lengths: ",plist_lengths)
  pli = fill(1,n)    # pli is a vector of indices for the array fields which is incremented lexigraplically
  max_pmap_list_length = 50000   # maximum number of elements in pmap_list
  count = 0
  result_dict = deepcopy(param_dict)
  for j = 1:n
    result_dict[array_field_list[j]] = param_dict[array_field_list[j]][pli[j]]
  end
  Base.push!( result, result_dict )
  i = 1
  while count < max_pmap_list_length
    i = 1
    while i <= n && pli[i] == plist_lengths[i]
      i += 1
    end
    if i <= n
      pli[i] += 1
      for j = 1:(i-1)
        pli[j] = 1
      end
    else
      break
    end
    result_dict = deepcopy(result_dict)
    for j = 1:n
      result_dict[array_field_list[j]] = param_dict[array_field_list[j]][pli[j]]
    end
    #println("pd: ",[ param_dict[array_field_list[j]][pli[j]] for j = 1:n]  )
    Base.push!( result, result_dict )
    count += 1
  end
  if count == max_pmap_list_length
    error("error:  count == max_pmap_list_length  in uni.jl.  Increase max_pmap_list_length.")
  end
  result
end
