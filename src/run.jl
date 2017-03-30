using TemporalEvolution

function run_trials( simname::AbstractString ) 
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  tr = temporal_result( N, num_attributes, num_subpops_list[1], ngens, mutation_stddev_list[1], num_emmigrants_list[1], 
      move_range, move_time_interval_list[1], opt_loss_cutoff, horiz_select_list[1], uniform_start )
  tr_list_run = TemporalEvolution.temporal_result_type[]
  trial=1
  for mutation_stddev in mutation_stddev_list
   for move_time_interval in move_time_interval_list
    for num_subpops in num_subpops_list
      for num_emmigrants in num_emmigrants_list
        for horiz_select in horiz_select_list
          tr = temporal_result( N, num_attributes, num_subpops, ngens, mutation_stddev, num_emmigrants, move_range, move_time_interval, 
              opt_loss_cutoff, horiz_select, uniform_start )
          Base.push!(tr_list_run, tr )
          #println("= = = = = = = =")
          #writerow(STDOUT,trial,tr)
        end
      end
    end
   end
  end
  println("===================================")
  tr_list_result = pmap(evolve, tr_list_run )
  trial = 1
  writeheader( stream, num_subpops_list, tr )
  writeheader( STDOUT, num_subpops_list, tr )
  for tr_result in tr_list_result
    writerow(stream,trial,tr_result)
    writerow(STDOUT,trial,tr_result)
    trial += 1
  end
end    

if length(ARGS) == 0
  simname = "configs/example2"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)
#println("simtype: ",simtype)
run_trials( simname )
