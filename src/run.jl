using TemporalEvolution

function run_trials( simname::AbstractString ) 
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  if !isdefined(:topology) && !isdefined(:topology_list)
    global topology="circular" 
  elseif isdefined(:topology) && !isdefined(:topology_list)
    println("isdef topology !isdef topology_list")
    global topology_list = [topology]
  end
  if !isdefined(:horiz_select)
    global horiz_select
    horiz_select = false
  end
  if !isdefined(:horiz_select_list)
    global horiz_select_list
    horiz_select_list = [false]
  end
  if !isdefined(:horiz_mutate)
    global horiz_mutate
    horiz_mutate = false
  end
  if !isdefined(:horiz_mutate_list)
    global horiz_mutate_list
    horiz_mutate_list = [false]
  end
  if !isdefined(:probHSelect)
    global probHSelect
    probHSelect = 1.0
  end
  if !isdefined(:probHSelect_list)
    global probHSelect_list
    probHSelect_list = [1.0]
  end
  println("simtype: ",simtype)
  println("topology_list: ",topology_list)
  #println("linear fitness: ",linear_fitness)
  int_burn_in = Int(round(burn_in*N))
  println("int_burn_in: ",int_burn_in)
  #horiz_mutate = false
  horiz_param_check( topology_list, num_subpops_list, num_emmigrants_list )
  tr = temporal_result( simtype, T, N, num_attributes_list[1], num_subpops_list[1], ngens, mutStddev_list[1], num_emmigrants_list[1], 
      move_range, move_time_interval_list[1], horiz_select, probHSelect_list[1], minFit, topology=topology_list[1], horiz_mutate=horiz_mutate_list[1],
      uniform_start=uniform_start, linear_fitness=linear_fitness, burn_in=burn_in, linfit_slope=linfit_slope )
  tr_list_run = TemporalEvolution.temporal_result_type[]
  trial=1
  for mutStddev in mutStddev_list
    for move_time_interval in move_time_interval_list
      for num_subpops in num_subpops_list
        for num_emmigrants in num_emmigrants_list
          for probHSelect in probHSelect_list
            for horiz_mutate in horiz_mutate_list
              for topology in topology_list
                for num_attributes in num_attributes_list
                  tr = temporal_result( simtype, T, N, num_attributes, num_subpops, ngens, mutStddev, num_emmigrants, move_range, move_time_interval, 
                      horiz_select, probHSelect, minFit, topology=topology, horiz_mutate=horiz_mutate,
                    uniform_start=uniform_start, linear_fitness=linear_fitness, burn_in=burn_in, linfit_slope=linfit_slope )
                  Base.push!(tr_list_run, tr )
                end
              end
            end
            #println("  length tr_list_run: ",length(tr_list_run))
          end
        end
      end
    end
  end
  println("===================================")
  tr_list_result = pmap(run_evolution, tr_list_run )
  #println("length tr_list_result: ",length(tr_list_result))
  #println("tr_list_result[1]: ",tr_list_result[1])
  #println("tr_list_result[1][1]: ",tr_list_result[1][1])
  #println("tr_list_result[1][2]: ",tr_list_result[1][2])
  trial = 1
  writeheader( stream, num_subpops_list, tr )
  writeheader( STDOUT, num_subpops_list, tr )
  for tr_result in tr_list_result
    writerow(stream,trial,tr_result)
    writerow(STDOUT,trial,tr_result)
    trial += 1
  end
  #=
  #println("trial: ",trial)
  close(stream)
  stream = open("$(simname)00.csv","w")
  writeheader( STDOUT, num_subpops_list, tr )
  #writeheader( stream, num_subpops_list, tr )
  for tr_res in tr_list_result
    for tr_result = tr_res[2]
      writerow(STDOUT,trial,tr_result)
      #writerow(stream,trial,tr_result)
      trial += 1
    end
  end
  close(stream)
  =#
end    

@everywhere function run_evolution( tr::temporal_result_type )
  if tr.simtype == 1
    return repeat_evolve_until_dead( tr )
  elseif tr.simtype == 2
    return repeat_evolve( tr )
  else
    error("illegal simtype in run.jl")
  end
end

if length(ARGS) == 0
  simname = "configs/example2"
else
  simname = ARGS[1]
  if length(ARGS) >= 2  # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
include("$(simname).jl")
println("simname: ",simname)
#println("simtype: ",simtype)
run_trials( simname )
