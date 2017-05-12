export print_temporal_result, writeheader, writerow, horiz_param_check
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run_spatial.jl configs/example1
=#
export  print_temporal_result, writeheader, writerow
#include("types.jl")
  
# Check topology and horizontal transfer settings for correctness
function horiz_param_check( topology_list::Vector{String}, num_subpops_list::Vector{Int64}, num_emmigrants_list::Vector{Int64} )
  #println("topology_list: ",topology_list)
  #println("num_subpops_list: ",num_subpops_list)
  #println("num_emmigrants_list: ",num_emmigrants_list)
  for topology in topology_list
    for num_subpops in num_subpops_list
      for ne in num_emmigrants_list
        #print("t: ",topology,"  nsbp: ",num_subpops,"  ne: ",ne,"  ")
        if topology == "none" || num_subpops == 1
          #println("OK")
        elseif ne > 0 && topology=="circular"
          #println("OK")
        elseif  ne > 0 && topology=="ring" || topology=="global"
          #println("OK")
        elseif  num_subpops >= 9 && ne > 0 && (topology=="vonneumann" || topology=="moore")
          #println("OK")
        elseif num_subpops > 1 && topology!="none" && ne > 0
          error("Warning!!! no horizontal transfer done with num_subpops=",num_subpops," and topology=",topology)
        end
      end
    end
  end
  println("horizontal transfer parameter check done.")
end

function print_temporal_result( tr::temporal_result_type )
  println("num_trials: ", tr.num_trials)
  println("N: ", tr.N)
  println("num_subpops: ", tr.num_subpops)
  println("num_emmigrants: ", tr.ne)
  println("num_attributes: ", tr.num_attributes)
  println("ngens: ", tr.ngens)
  println("burn_in: ", tr.burn_in)
  println("uniform_start: ", tr.uniform_start)
  println("horiz_select: ", tr.horiz_select)
  println("probHSelect: ", tr.probHSelect)
  println("mutStddev: ", tr.mutStddev)
  println("ideal_init: ", tr.ideal_init)
  println("move_range: ", tr.move_range)
  println("move_time_interval: ", tr.move_time_interval)
  println("linear_fitness: ",tr.linear_fitness)
  println("linfit_slope: ",tr.linfit_slope)
  println("topology: ",tr.topology)
  println("minFit: ",tr.minFit)
  println("mean_fraction_subpops_below_cutoff: ", tr.mean_fraction_subpops_below_minFit)
  println("fraction_gens_with_all_subpops_below_minFit: ", tr.fraction_gens_with_all_subpops_below_minFit)
  println("fitness_mean: ", tr.fitness_mean)
  println("fitness_variance: ", tr.fitness_variance)
  println("attiribute_variance: ", tr.attribute_variance)
  println("generational_lifetime: ", tr.generational_lifetime)
  println("move_update_lifetime: ", tr.move_update_lifetime)
  println("gen_limit_reached_count: ", tr.gen_limit_reached_count)
end

function writeheader( stream::IO, num_subpops_list::Vector{Int64}, tr::temporal_result_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# num_trials=$(tr.num_trials)",
    "# N=$(tr.N)",
    "# num_subpops_list=$(num_subpops_list)",
    #"# num_attributes=$(tr.num_attributes)",
    #"# horiz_select=$(tr.horiz_select)",
    #"# num_emmigrants=$(tr.ne)",
    "# ngens=$(tr.ngens)",
    "# burn_in=$(tr.burn_in)",
    #"# mutStddev=$(tr.mutStddev)",
    "# uniform_start=$(tr.uniform_start)",
    "# minFit=$(tr.minFit)",
    "# linear_fitness=$(tr.linear_fitness)",
    #"# linfit_slope=$(tr.linfit_slope)",
    #"# topology=$(tr.topology)"
  ]
  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "num_subpops",
    "subpop_size",
    "num_emigrants",
    "mutStddev",
    "num_attributes",
    "move_range",
    "move_time_interval",
    "horiz_select",
    "topology",
    "linear_fitness_slope"
    ]
  if tr.simtype == 2
    heads2 = [
    "mean_fitness",
    "stddev_fitness",
    "stddev_attributes",
    "mean_fraction_subpops_below_minFit",
    "fraction_gens_with_all_subpops_below_minFit"
    ]
  elseif tr.simtype == 1
    heads2 = [
    "generational_lifetime",
    "move_update_lifetime",
    "gen_limit_reached_count"
    ]
  end
  append!(heads,heads2)
  write(stream,join(heads,","),"\n")
end
    
function writerow( stream::IO, trial::Int64, tr::temporal_result_type )
  line = Any[
      tr.num_subpops,
      Int(floor(tr.N/tr.num_subpops)),
      tr.ne,
      tr.mutStddev,
      tr.num_attributes,
      tr.move_range,
      tr.move_time_interval,
      tr.horiz_select,
      tr.topology,
      tr.linfit_slope
  ]
  if tr.simtype==2
    line2 = Any[
      tr.fitness_mean,
      sqrt(tr.fitness_variance),
      sqrt(tr.attribute_variance),
      tr.mean_fraction_subpops_below_minFit,
      tr.fraction_gens_with_all_subpops_below_minFit,
    ]
  elseif tr.simtype == 1
    line2 = Any[
      tr.generational_lifetime,
      tr.move_update_lifetime,
      tr.gen_limit_reached_count
    ]
  end  
  append!(line,line2)
  write(stream,join(line,","),"\n")
end

