export print_temporal_result, writeheader, writerow
#=
Recommended command line to run:
>  julia -L TemporalEvolution.jl run_spatial.jl configs/example1
=#
export  print_temporal_result, writeheader, writerow
#include("types.jl")
  

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
  println("mutation_stddev: ", tr.mutation_stddev)
  println("ideal_init: ", tr.ideal_init)
  println("move_range: ", tr.move_range)
  println("move_time_interval: ", tr.move_time_interval)
  println("opt_loss_cutoff: ", tr.opt_loss_cutoff)
  println("linear_fitness: ",tr.linear_fitness)
  println("topology: ",tr.topology)
  println("min_fit: ",tr.min_fit)
  println("mean_fraction_subpops_below_cutoff: ", tr.mean_fraction_subpops_below_cutoff)
  println("fitness_mean: ", tr.fitness_mean)
  println("fitness_variance: ", tr.fitness_variance)
  println("attiribute_variance: ", tr.attribute_variance)
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
    #"# mutation_stddev=$(tr.mutation_stddev)",
    "# opt_loss_cutoff=$(tr.opt_loss_cutoff)",
    "# uniform_start=$(tr.uniform_start)",
    "# min_fit=$(tr.min_fit)",
    "# linear_fitness=$(tr.linear_fitness)",
    "# topology=$(tr.topology)"
  ]
  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "num_subpops",
    "subpop_size",
    "num_emigrants",
    "mutation stddev",
    "num_attributes",
    "move_range",
    "move_time_interval",
    "horiz_select",
    "mean_fitness",
    "stddev_fitness",
    "stddev_attributes",
    "fract_below_cutoff"
  ]
  write(stream,join(heads,","),"\n")
end
    
function writerow( stream::IO, trial::Int64, tr::temporal_result_type )
  line = Any[
          tr.num_subpops,
          Int(floor(tr.N/tr.num_subpops)),
          tr.ne,
          tr.mutation_stddev,
          tr.num_attributes,
          tr.move_range,
          tr.move_time_interval,
          tr.horiz_select,
          tr.fitness_mean,
          sqrt(tr.fitness_variance),
          sqrt(tr.attribute_variance),
          tr.mean_fraction_subpops_below_cutoff
  ]
  write(stream,join(line,","),"\n")
end

