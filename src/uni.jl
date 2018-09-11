# Legal fields for the dictionaries that are used to parameters, results, and individuals.
# Suggested test sequence:
# julia> paramd = init_dictionary( param_fields )   # builds dictionary paramd
# julia> read_parameter_file("param",paramd)        # adds parameters from param file into dictionary
# julia> build_paramd_list( paramd )                  # builds list of dictionaries for each trial

export temporal_param_fields, init_dictionary, read_parameter_file, print_meta_pop, build_paramd_list, 
      print_meta_pop_attributes, count_alive, subpops_alive, opt_gained_lost_counts

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

function print_meta_pop( paramd::param_type, meta_pop::PopList, vt::Dict{Int64,temporal_variant_type} )
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

# counts number of alive individuals in subpop
# An individual is "alive" if it has fitness above minFit
# An individual is "dead" if it has fitness == minFit
function count_alive( subpop::Population, paramd::param_type, vt::Dict{Int64,temporal_variant_type})
  #print("count_alive: ")
  count = 0
  for x in subpop
    #print(" vt[",x,"] fitness: ",vt[x].fitness)
    if vt[x].fitness > paramd[:minFit]
      count += 1
    end
  end
  #println("  count: ",count)
  return count
end

# Returns a Bool array over subpops where true means alive 
function subpops_alive( meta_pop::PopList, paramd::param_type, vt::Dict{Int64,temporal_variant_type})
  result = fill(false, length(meta_pop) )
  i = 1
  for s in meta_pop
    if count_alive( s, paramd, vt ) > 0
      result[i] = true
    end
    i += 1
  end
  #println("subpop_alive: ",result)
  result
end

# Each of the parameters is a Boolean vector over subpopulations
# Returns the count of subpops for the following situations
# subpop gains optimum after move_optima and propsel
# subpop loses optimum after move_optima and propsel
# subpop gains optimum after horiz
# subpop loses optimum after horiz
# Note that "gains optimum" is equivalent going from dead to alive.
# Note that "loses optimum" is equivalent going from alive to dead.
function opt_gained_lost_counts( prev_gen::Vector{Bool}, after_propsel::Vector{Bool}, after_horiz::Vector{Bool} )
  (propsel_loss,propsel_gain,horiz_loss,horiz_gain) = (0,0,0,0)
  num_subpops = length(prev_gen)
  for i = 1:num_subpops
    propsel__loss = (prev_gen[i] && !after_propsel[i]) ? 1 : 0
    propsel__gain = (!prev_gen[i] && after_propsel[i]) ? 1 : 0
    horiz__loss = (after_propsel[i] && !after_horiz[i]) ? 1 : 0
    horiz__gain = (!after_propsel[i] && after_horiz[i]) ? 1 : 0
    #println("i: ",i," ",(propsel__loss,propsel__gain,horiz__loss,horiz__gain))
    (propsel_loss,propsel_gain,horiz_loss,horiz_gain) = 
        (propsel_loss,propsel_gain,horiz_loss,horiz_gain) .+  (propsel__loss,propsel__gain,horiz__loss,horiz__gain) 
  end
  (propsel_loss,propsel_gain,horiz_loss,horiz_gain)
end
  

function print_meta_pop_attributes( paramd::param_type, meta_pop::PopList, vt::Dict{Int64,temporal_variant_type} )
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
    if param_dict[key] == :null && key != :num_fit_locations
      println("warning:  parameter: ",key," not initialized")
    end
  end
  param_dict
end  

function check_parameters!( paramd::Dict{Symbol,Any} )
  if paramd[:num_subpops] == "sqrt"
    num_subpops = Int(floor(N/sqrt(N )))
    #subpop_size = Int(floor(sqrt( N )))
    subpop_size = Int(floor(N/num_subpops))
    paramd[:N] = num_subpops*subpop_size
    println("N: ",paramd[:N]N,"subpop_size: ",subpop_size,"  num_subpops: ",num_subpops)
  end
end

# builds a list of dictionaries that can be used to run a sequence of trials with different
#   parameter settings using pmap()  (parallel map).
function build_paramd_list(  paramd::Dict{Symbol,Any} )
  #println("build_paramd_list: N: ",paramd[:N])
  array_field_list = Symbol[]
  for k in keys(paramd)
    if typeof(paramd[k]) <: Array
      Base.push!( array_field_list, k )
    end
  end
  sort!(array_field_list,rev=true)  # use reverse sorted order to be compatible with spatial
  paramd_list = []
  n = length(array_field_list)
  plist_lengths = [ length( paramd[k] ) for k in array_field_list ]
  #println("plist_lengths: ",plist_lengths)
  pli = fill(1,n)    # pli is a vector of indices for the array fields which is incremented lexigraplically
  max_paramd_list_length = 50000   # maximum number of elements in paramd_list
  count = 0
  paramd_list_dict = deepcopy(paramd)
  #i = 1
  while count < max_paramd_list_length
    paramd_list_dict = deepcopy(paramd_list_dict)
    for j = 1:n
      paramd_list_dict[array_field_list[j]] = paramd[array_field_list[j]][pli[j]]
    end
    # Handle the special situation where paramd[:num_subpops] == "sqrt"
    if paramd[:num_subpops] == "sqrt"
      paramd_list_dict[:num_subpops] = Int(floor(paramd_list_dict[:N]/sqrt( paramd_list_dict[:N] )))
      paramd_list_dict[:subpop_size] = Int(floor( paramd_list_dict[:N]/paramd_list_dict[:num_subpops] ))
      paramd_list_dict[:N] = paramd_list_dict[:num_subpops]*paramd_list_dict[:subpop_size]
    else
      paramd_list_dict[:subpop_size] = Int(floor( paramd_list_dict[:N]/paramd_list_dict[:num_subpops] ))
    end
    # Handle the special spatial computation of paramd[:num_fit_locations]
    if paramd[:simtype] == 4  # Spatial
      if typeof(paramd[:num_subpops]) <: Array
        paramd_list_dict[:num_fit_locations] = paramd_list_dict[:use_fit_locations] ? maximum(paramd[:num_subpops]) : paramd_list_dict[:num_subpops]
      else
        paramd_list_dict[:num_fit_locations] = paramd_list_dict[:num_subpops]
      end
    end
    #print_dict("uni: paramd_list_dict:",paramd_list_dict)
    Base.push!( paramd_list, paramd_list_dict )
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
    count += 1
  end
  if count == max_paramd_list_length
    error("error:  count == max_paramd_list_length  in uni.jl.  Increase max_paramd_list_length.")
  end
  paramd_list
end
