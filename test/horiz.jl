include("../src/TemporalEvolution.jl")
#include("../src/ev_init.jl")
function test_new_emmigrants()
  #srand(1)
  tr = ev_init()
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  ideal = fill( tr.ideal_init, tr.num_attributes )
  vt = Dict{Int64,variant_type}()
  if tr.uniform_start
    meta_pop = [  fill(1,subpop_size)  for j = 1:tr.num_subpops ]
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tr.num_subpops ]
    for i = 1:tr.N
      attr = rand(tr.num_attributes)
      vt[i] = variant_type( fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness ), attr )
    end
  end
  #println("vt: ",vt)
  neighbor_list = [ (i % tr.num_subpops) + 1 for i = 1:num_subpops ]
  id = Int64[0]
  for i = 1:tr.N
    attr = rand(tr.num_attributes)
    vt[i] = variant_type( TemporalEvolution.fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness ), attr )
  end
  id[1] += tr.N
  new_emmigrants = new_emmigrants_funct( meta_pop, tr, vt, neighbor_list, ideal, id, emmigrant_select=false )
  meta_pop= add_emmigrants( meta_pop, tr, vt, new_emmigrants, neg_select=horiz_select )
  println("meta_pop: ",meta_pop)
end
#test_new_emmigrants()

function test_horiz_transfer_by_fitness()
  srand(1)
  tr = ev_init()
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  ideal = fill( tr.ideal_init, tr.num_attributes )
  vt = Dict{Int64,variant_type}()
  if tr.uniform_start
    meta_pop = [  fill(1,subpop_size)  for j = 1:tr.num_subpops ]
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tr.num_subpops ]
    for i = 1:tr.N
      attr = rand(tr.num_attributes)
      vt[i] = variant_type( fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness ), attr )
    end
  end
  #println("vt: ",vt)
  neighbor_list = [ (i % tr.num_subpops) + 1 for i = 1:num_subpops ]
  id = Int64[0]
  for i = 1:tr.N
    attr = rand(tr.num_attributes)
    vt[i] = variant_type( TemporalEvolution.fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness ), attr )
  end
  (mmeans, vvars) = means( meta_pop, vt )
  id[1] += tr.N
  horiz_transfer_by_fitness!( meta_pop, tr, vt, ideal, mmeans, id, topology="vonneumann", neg_select=horiz_select, emmigrant_select=horiz_select )
  println("meta_pop: ",meta_pop)
end

test_horiz_transfer_by_fitness()
