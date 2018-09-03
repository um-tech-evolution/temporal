export print_temporal_result, writeheader, writerow, horiz_param_check
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run_spatial.jl configs/example1
=#
export  print_temporal_result, writeheader, writerow, print_dict
#include("types.jl")
  
# TODO:  Automate this function instead of hard coding  the fields
function print_temporal_params( tp::param_type )
  println("simtype: ", tp[:simtype])
  println("num_trials: ", tp[:num_trials])
  println("N: ", tp[:N])
  println("num_subpops: ", tp[:num_subpops])
  println("num_emigrants: ", tp[:num_emigrants])
  println("num_attributes: ", tp[:num_attributes])
  println("ngens: ", tp[:ngens])
  println("burn_in: ", tp[:burn_in])
  println("uniform_start: ", tp[:uniform_start])
  println("horiz_select: ", tp[:horiz_select])
  println("probHSelect: ", tp[:probHSelect])
  println("mutStddev: ", tp[:mutStddev])
  println("ideal_init: ", tp[:ideal_init])
  println("move_range: ", tp[:move_range])
  println("move_time_interval: ", tp[:move_time_interval])
  println("linfit_slope: ",tp[:linfit_slope])
  println("horiz_mutate: ",tp[:horiz_mutate])
  println("topology: ",tp[:topology])
  println("minFit: ",tp[:minFit])
end 

# TODO:  Automate this function instead of hard coding  the fields
function print_temporal_results( tr::result_type )
  println("mean_fraction_subpops_below_cutoff: ", tr[:mean_fraction_subpops_below_minFit])
  println("fraction_gens_with_all_subpops_below_minFit: ", tr[:fraction_gens_with_all_subpops_below_minFit])
  println("fitness_mean: ", tr[:fitness_mean])
  println("fitness_variance: ", tr[:fitness_variance])
  println("attiribute_variance: ", tr[:attribute_variance])
  println("generational_lifetime: ", tr[:generational_lifetime])
  println("move_update_lifetime: ", tr[:move_update_lifetime])
  println("gen_limit_reached_count: ", tr[:gen_limit_reached_count])
  println("innovations_per_gen_trial: ", tr[:innovations_per_gen_trial])
  println("half_innovations_per_gen_trial: ", tr[:half_innovations_per_gen_trial])
  println("deleterious_mutations_per_gen_trial: ", tr[:innovations_per_gen_trial])
  println("half_deleterious_mutations_per_gen_trial: ", tr[:half_innovations_per_gen_trial])
end

function print_dict( name::String, dict::Dict{Symbol,Any} )
  println(name)
  for key in keys(dict)
    println(key,": ",dict[key])
  end
end

@doc """ function writeheader() 
  paramd is the dictionary of parameters as read in from the configuration file.
  These parameters are written one per line preceded by the comment character hash.
  resultd is the dictionary that includes those parameters that change from trial to trial, plus the results
     from each trial.  Each key of this dictionary corresponds to a column of the output CSV file.  This
     function writes the headers for the columns of the CSV file.
"""
function writeheader( stream::IO, paramd::param_type, resultd::result_type )
  seed = try seed = Main.seed catch -1 end  # assign local variable seed to the Main module seed if defined, otherwise to -1 
  param_fields = paramd[:simtype] == 4 ? spatial_param_fields : temporal_param_fields
  param_strings = String["# $(string(Dates.today()))"]
  Base.push!(param_strings,(seed!==-1 ? "# random number seed: $(seed)": "# no random number seed defined"))
  for k in param_fields
    if paramd[k] != :null
      Base.push!( param_strings, "# $(String(k))=$(paramd[k])")
    end
  end
  write(stream,join(param_strings,"\n"),"\n")   
  heads = String[]
  for k in param_fields
    if paramd[k] != :null && typeof(paramd[k]) <: Array 
      #println("writeheader: array field: ",k,paramd[k])
      Base.push!(heads,String(k))
    end
  end 
  if  !(typeof(paramd[:num_subpops]) <: Array) 
    Base.push!(heads,String(:num_subpops))   # always include the :num_subpops field
  end
  Base.push!(heads,String(:subpop_size))   # always include the :subpop_size field
  for k in keys(resultd)
    Base.push!(heads,String(k))
  end
  write(stream,join(heads,","),"\n")
end
    
# paramd  is the parameter dictionary read in from the parameter file which includes list (or array) values
# tp  is the parameter dictionary corresponding to a particular trial.  All values are fully specified.
function writerow( stream::IO, paramd::param_type, tp::param_type, resultd::result_type )
  #println("writerow paramd[:num_subpops]: ",paramd[:num_subpops],"  topology: ",paramd[:topology])
  param_fields = paramd[:simtype] == 4 ? spatial_param_fields : temporal_param_fields
  line = Any[]
  #for k in keys(paramd)
  for k in param_fields
    if paramd[k] != :null && typeof(paramd[k]) <: Array
      #println("writerow: array field: ",k,paramd[k])
      Base.push!(line,tp[k])
    end
  end 
  if  !(typeof(paramd[:num_subpops]) <: Array) 
    Base.push!(line,tp[:num_subpops])   # always include the :num_subpops field
  end
  Base.push!(line,tp[:subpop_size])   # always include the :subpop_size field
  for k in keys(resultd)
    Base.push!(line,resultd[k])
  end 
  write(stream,join(line,","),"\n")
end

