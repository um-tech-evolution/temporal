@everywhere include("TemporalEvolution.jl")
#include("uni.jl")
#include("types.jl")
export seed

#= Dictionaries:
In this program, dictionaries are used to store composite data structures.
The relevant dictionaries and lists of dictionaries are:
 *  paramd:  read from the configuration file.  
             Dictionary keys are parameters and dictionary values are parameter values
                (which may be lists if there is iteration over multiple paramter values)
 *  resultd: Dictionary keys are the results from each trial.  Values are null.
 *  paramd_list:  List of parameter dictionaries, one for each trial, with the same fields as paramd.
                  Now parameters have the specified values used for the corresponding trial.
 *  resultd_list:  List of result dictionaries, one for each trial, with the same fields as resultd.
                  The values are the computed results for the corresponding trial.
There is one special situation, namely when the parameter :num_subpops has the value "sqrt" rather
than a numerical value in the configuration file and in paramd.  In this case, the paramd_list value
of :num_subpops is set to Int(floor(sqrt(N))) where N is the meta_pop size for the corresponding trial.  
Then the :N value is rounded down to a perfect square, and :subpop_size is set to :num_subpops.

=#

@everywhere function run_evolution( paramd::TemporalEvolution.param_type, resultd::TemporalEvolution.result_type )
  #println("run_evolution: simtype: ",paramd[:simtype])
  if paramd[:simtype] == 1
    return repeat_evolve_until_dead( paramd, deepcopy(resultd) )
  elseif paramd[:simtype] == 2
    #println("run_evolution repeat_evolve")
    return repeat_evolve( paramd, deepcopy(resultd) )
  elseif paramd[:simtype] == 4
    #println("run_evolution: num_fit_locations: ",paramd[:num_fit_locations])
    return repeat_spatial( paramd, deepcopy(resultd) )
  else
    error("illegal simtype in run.jl")
  end
end

function run_trials( paramd::TemporalEvolution.param_type, resultd::TemporalEvolution.result_type )
  stream = open("$(paramd[:simname]).csv","w")
  paramd_list = build_paramd_list( paramd )
  resultd_list = pmap(x->run_evolution(x,resultd), paramd_list )     #  comment out for debugging: error messages are simpler
  #resultd_list = map(x->run_evolution(x,resultd), paramd_list )     # TODO uncomment for debugging: error messages are simpler
  r = resultd_list[1]
  if paramd[:simtype] == 2
    #println("fitness_mean: ",r[:fitness_mean],"  attr variance: ",r[:attribute_variance])
  elseif paramd[:simtype] == 1
    #println("generational_lifetime: ",r[:generational_lifetime])
  elseif paramd[:simtype] == 4
    #println("fitness_mean: ",r[:fitness_mean],"  attr variance: ",r[:attribute_variance])
  end
  #print_dict("Dictionary paramd",paramd)
  #print_dict("Dictionary resultd",resultd)
  #print_dict("Dictionary resultd_list[1]: ",resultd_list[1])
  writeheader( STDOUT, paramd, resultd )
  writeheader( stream, paramd, resultd )
  for i = 1:length(paramd_list)
    writerow( STDOUT, paramd, paramd_list[i], resultd_list[i] )
    writerow( stream, paramd, paramd_list[i], resultd_list[i] )
  end
end

global seed
global param_fields_list, result_fields_list
if length(ARGS) == 0
  simname = "examples/example1"
else
  simname = ARGS[1]
  if length(ARGS) >= 2  # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
println("simname: ",simname)
include("$(simname).jl")

if simtype == 2
  println("simtype == 2 which is temporal reporting equilibrium stats over generations.")
  param_fields_list = TemporalEvolution.temporal_param_fields   
  result_fields_list = TemporalEvolution.temporal_result_fields
elseif simtype == 1
  println("simtype == 1 which is temporal reporting number of generations until ideal is lost")
  param_fields_list =  TemporalEvolution.temporal_param_fields
  result_fields_list =  TemporalEvolution.subpop_alive_result_fields
elseif simtype == 4
  println("simtype == 4 which is spatial reporting equilibrium stats over generations.")
  param_fields_list = TemporalEvolution.spatial_param_fields 
  result_fields_list = TemporalEvolution.spatial_result_fields 
end
paramd =init_dictionary( param_fields_list )
paramd[:simname] = simname
paramd = read_parameter_file( "$(simname).jl", paramd )
resultd = init_dictionary( result_fields_list )
run_trials( paramd, resultd )
