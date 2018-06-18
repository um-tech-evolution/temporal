@everywhere include("TemporalEvolution.jl")
include("uni.jl")
include("types.jl")


@everywhere function run_evolution( paramd::TemporalEvolution.param_type, resultd::TemporalEvolution.result_type )
  if paramd[:simtype] == 1
    return repeat_evolve_until_dead( paramd, deepcopy(resultd) )
  elseif paramd[:simtype] == 2
    return repeat_evolve( paramd, deepcopy(resultd) )
  else
    error("illegal simtype in run.jl")
  end
end

function run_trials( paramd::TemporalEvolution.param_type, resultd::TemporalEvolution.result_type )
  stream = open("$(paramd[:simname]).csv","w")
  pmap_list = build_pmap_list( paramd )
  resultd_list = pmap(x->run_evolution(x,resultd), pmap_list )
  r = resultd_list[1]
  if paramd[:simtype] == 2
    #println("fitness_mean: ",r[:fitness_mean],"  attr variance: ",r[:attribute_variance])
  else
    #println("generational_lifetime: ",r[:generational_lifetime])
  end
  writeheader( STDOUT, paramd, resultd )
  writeheader( stream, paramd, resultd )
  for i = 1:length(pmap_list)
    writerow( STDOUT, paramd, pmap_list[i], resultd_list[i] )
    writerow( stream, paramd, pmap_list[i], resultd_list[i] )
  end
end

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
paramd = init_dictionary( TemporalEvolution.temporal_param_fields )
paramd[:simname] = simname
paramd = read_parameter_file( "$(simname).jl", paramd )
println("simtype: ",paramd[:simtype])
if paramd[:simtype] == 2
  resultd = init_dictionary( TemporalEvolution.temporal_result_fields )
elseif paramd[:simtype] == 1
  resultd = init_dictionary( TemporalEvolution.subpop_alive_result_fields )
end
run_trials( paramd, resultd )
