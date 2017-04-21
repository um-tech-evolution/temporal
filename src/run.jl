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
  println("topology_list: ",topology_list)
  #println("linear fitness: ",linear_fitness)
  #println("burn_in: ",burn_in)
  tr = temporal_result( T, N, num_attributes, num_subpops_list[1], ngens, mutation_stddev_list[1], num_emmigrants_list[1], 
      move_range, move_time_interval_list[1], horiz_select_list[1], min_fit, topology=topology_list[1],
      uniform_start=uniform_start, linear_fitness=linear_fitness, burn_in=burn_in, linfit_slope=linfit_slope_list[1] )
  tr_list_run = TemporalEvolution.temporal_result_type[]
  trial=1
  for mutation_stddev in mutation_stddev_list
    for move_time_interval in move_time_interval_list
      for num_subpops in num_subpops_list
        for num_emmigrants in num_emmigrants_list
          for horiz_select in horiz_select_list
            for topology in topology_list
              for linfit_slope in linfit_slope_list
                tr = temporal_result( T, N, num_attributes, num_subpops, ngens, mutation_stddev, num_emmigrants, move_range, move_time_interval, 
                  horiz_select, min_fit, topology=topology,
                  uniform_start=uniform_start, linear_fitness=linear_fitness, burn_in=burn_in, linfit_slope=linfit_slope )
                Base.push!(tr_list_run, tr )
              end
            end
            #println("  length tr_list_run: ",length(tr_list_run))
          end
        end
      end
    end
  end
  println("===================================")
  tr_list_result = pmap(repeat_evolve, tr_list_run )
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
  writeheader( stream, num_subpops_list, tr )
  for tr_res in tr_list_result
    for tr_result = tr_res[2]
      writerow(STDOUT,trial,tr_result)
      writerow(stream,trial,tr_result)
      trial += 1
    end
  end
  close(stream)
  =#
end    

if length(ARGS) == 0
  simname = "configs/example2"
else
  simname = ARGS[1]
end
#srand(1)
include("$(simname).jl")
println("simname: ",simname)
#println("simtype: ",simtype)
run_trials( simname )
