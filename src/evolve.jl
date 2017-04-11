export temporal_result, repeat_evolve, evolve, mutate_meta_pop!, means
    
#include("types.jl")
#include("propsel.jl")
#include("fitness.jl") 

function temporal_result( num_trials::Int64, N::Int64, num_attributes::Int64, num_subpops::Int64, ngens::Int64, mutation_stddev::Float64, num_emmigrants::Int64,
    move_range::Float64, move_time_interval::Int64, opt_loss_cutoff::Float64, horiz_select::Bool=false, min_fit::Float64=0.0; topology::String="circular", 
    uniform_start::Bool=false, linear_fitness::Bool=false, burn_in::Float64=1.0 )
  ideal_init = 0.5
  return temporal_result_type( num_trials, N, num_subpops, num_emmigrants, num_attributes, ngens, burn_in, uniform_start, horiz_select, mutation_stddev,
      ideal_init, move_range, move_time_interval, opt_loss_cutoff, min_fit, linear_fitness, topology, 0.0, 0.0, 0.0, 0.0 )
end

@doc """ function repeat_evolve( )
Runs evolve()  tr.num_subpops times, and averages the results.
"""
function repeat_evolve( tr::temporal_result_type )
  println("rep_evolve: num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutation stddev: ",tr.mutation_stddev)
  if tr.num_trials == 1
    return evolve( tr )
  end
  tr_list = temporal_result_type[]
  sum_mean = 0.0
  sum_vars = 0.0
  sum_attr_vars = 0.0
  sum_count_below_cutoff = 0

  for t = 1:tr.num_trials
    tr = evolve( tr ) 
    Base.push!( tr_list, deepcopy(tr) )
    sum_count_below_cutoff += tr.mean_fraction_subpops_below_cutoff
    #println("t: ",t," tr.fitness_variance: ",tr.fitness_variance)
    sum_mean += tr.fitness_mean
    sum_vars += tr.fitness_variance
    sum_attr_vars += tr.attribute_variance
  end
  tr.mean_fraction_subpops_below_cutoff = sum_count_below_cutoff/tr.num_trials
  tr.fitness_mean = sum_mean/tr.num_trials
  tr.fitness_variance = sum_vars/tr.num_trials
  #println("  tr.fitness_variance: ",tr.fitness_variance)
  #println("  tr.num_trials: ",tr.num_trials)
  tr.attribute_variance = sum_attr_vars/tr.num_trials
  #return tr, tr_list
  return tr
end

@doc """ function evolve( )
The main simulation function for temporal.
See types.jl for the definition of temporal_result_type, and for the definition of the fields of this type.
"""
function evolve( tr::temporal_result_type )
  #println("num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutation stddev: ",tr.mutation_stddev,"  topology: ",tr.topology)
  int_burn_in = Int(round(tr.burn_in*tr.N))
  id = [0]
  mmeans = zeros(tr.num_subpops)
  vvars = zeros(tr.num_subpops)
  att_vars = zeros(tr.num_subpops)
  cumm_means = zeros(tr.num_subpops)
  cumm_vars = zeros(tr.num_subpops)
  cumm_attr_vars = zeros(tr.num_subpops)
  count_below_cutoff = 0
  cumm_count_below_cutoff = 0
  ideal = fill( tr.ideal_init, tr.num_attributes )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  if subpop_size*tr.num_subpops != tr.N
    error("N must equal subpop_size*tr.num_subpops")
  end
  vt = Dict{Int64,variant_type}()
  if tr.uniform_start
    meta_pop = [  fill(1,subpop_size)  for j = 1:tr.num_subpops ]
    attr = ideal
    vt[1] = variant_type( fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness ), attr )
    id[1] += 1
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tr.num_subpops ]
    for i = 1:tr.N
      attr = rand(tr.num_attributes)
      vt[i] = variant_type( fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness ), attr )
    end
    id[1] += tr.N
  end
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  for g = 2:(tr.ngens+int_burn_in)
    if g > int_burn_in && g % tr.move_time_interval == 0
      #println("g: ",g,"  mmeans: ",mmeans)
      #print("   fstdev: ",sqrt(vvars))
      #print("   astdev: ",sqrt(att_vars))
      #println("  count below cutoff: ",count_pops_below_fit_cutoff( mmeans, tr.opt_loss_cutoff ))
      #println("  count below cutoff: ",count_below_cutoff,"   cumm_count_below_cutoff: ",cumm_count_below_cutoff)
      move_optima( ideal, tr.move_range )
    end
    mutate_meta_pop!( meta_pop, vt, ideal, id, tr )  # will also re-evaluate fitness
    for  j = 1:tr.num_subpops
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )
    end
    mmeans, vvars = means( meta_pop, vt )
    if tr.num_subpops >= 9 && tr.ne > 0 && tr.topology=="circular"
      horiz_transfer_circular!( meta_pop, tr, vt, ideal, id, g, neg_select=tr.horiz_select, emmigrant_select=tr.horiz_select )
    elseif  tr.num_subpops >= 9 && tr.ne > 0 && tr.topology!="circular"
      horiz_transfer_by_fitness!( meta_pop, tr, vt, ideal, mmeans, id, neg_select=tr.horiz_select, topology=tr.topology, emmigrant_select=tr.horiz_select )
    end
    if g > int_burn_in  # data collection
      #(mmeans, vvars) = means( meta_pop, vt )
      cumm_means += mmeans
      cumm_vars += vvars
      #print("g: ",g,"  mmeans: ",mmeans)
      #print("   fstdev: ",sqrt(vvars))
      att_vars = attr_vars( meta_pop, vt )
      cumm_attr_vars += att_vars
      #println("   astdev: ",sqrt(att_vars))
      count_below_cutoff = count_pops_below_fit_cutoff( mmeans, tr.opt_loss_cutoff )
      cumm_count_below_cutoff += count_below_cutoff
      #println("  count below cutoff: ",count_below_cutoff,"   cumm_count_below_cutoff: ",cumm_count_below_cutoff)
    end
  end
  #tr.mean_fraction_subpops_below_cutoff = cumm_count_below_cutoff/tr.num_subpops/tr.ngens  # commented out 4/4/17
  tr.mean_fraction_subpops_below_cutoff = cumm_count_below_cutoff/tr.ngens  # added  4/4/17
  tr.fitness_mean = mean(cumm_means/tr.ngens)
  #println("cumm_means/tr.ngens: ",cumm_means/tr.ngens,"  mean: ",tr.fitness_mean)
  #println("cumm_count_below_cutoff: ",cumm_count_below_cutoff)
  #println("cumm_count_below_cutoff/num_subpops: ",cumm_count_below_cutoff/tr.num_subpops)
  tr.fitness_variance = mean(cumm_vars/tr.ngens)
  tr.attribute_variance = mean(cumm_attr_vars/tr.ngens)
  return tr
end

function mutate_meta_pop!( meta_pop::PopList, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64}, tr::temporal_result_type  )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  #println("num_subpops: ",num_subpops,"  subpop_size: ",subpop_size)
  #v_lists = [ [ vt[meta_pop[j][i]] for i = 1:subpop_size ] for j = 1:tr.num_subpops]
  v_lists = [ [ vt[meta_pop[j][i]] for i = 1:subpop_size ] for j = 1:num_subpops ]
  for j = 1:num_subpops
    #v_list = [ vt[meta_pop[j][i]] for i = 1:subpop_size ]
    for i = 1:subpop_size
      #println("B j: ",j,"  i: ",i," meta_pop[j][i]: ",meta_pop[j][i],"  v_lists[j][i]: ",v_lists[j][i])
      meta_pop[j][i] = id[1]
      id[1] += 1
      vt[meta_pop[j][i]] = mutate_variant( v_lists[j][i], tr.mutation_stddev )
      vt[meta_pop[j][i]].fitness = fitness( vt[meta_pop[j][i]].attributes, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness )
      #println("A j: ",j,"  i: ",i," meta_pop[j][i]: ",meta_pop[j][i],"  vt[meta_pop[j][i]]: ",vt[meta_pop[j][i]])
    end
  end
end


function mutate_variant( v::variant_type, mutation_stddev::Float64 )
  new_v = deepcopy(v)
  new_v.attributes = mutate_attributes( new_v.attributes, mutation_stddev )
  new_v
end

function mutate_attributes( attributes::Vector{Float64}, mutation_stddev::Float64 )
  #println("mutate attributes  mutation_stddev: ",mutation_stddev)
  for i = 1:length(attributes)
    #println("B attributes[",i,"]: ",attributes[i])
    attributes[i] += +mutation_stddev*randn()
    while attributes[i] < 0
        attributes[i] += 1.0
        #println("wrapped up: ",attributes[i])
    end
    while attributes[i] > 1.0
        attributes[i] -= 1.0
        #println("wrapped down: ",attributes[i])
    end
    if attributes[i] < 0.0 || attributes[i] > 1.0
      error("attribute not wrapped")
    end
    attributes[i] = min(1.0,max(0.0,attributes[i]))
    #println("A attributes[",i,"]: ",attributes[i])
  end
  #println("attributes: ",attributes)
  return attributes
end

@doc """ move_optima()
  Add a float between -move_range and +move_range to each ideal value.
  If discrete_move==true, then the float is randomly chosen in this range.
  If discrete_move==false, then the float is randomly chosen to be either +move_range or -move_range.
"""
function move_optima( ideal::Vector{Float64}, move_range::Float64; discrete_move::Bool=true )
  num_attributes = length(ideal)
  for k = 1:num_attributes
    if discrete_move   # Move by move_range or -move_range  
      move = rand() > 0.5 ? move_range : -move_range
    else  # move by a randomly chosen amount between -move_range and +move_range
      move = 2.0*rand()*move_range - move_range
    end
    ideal[k] += move
    while ideal[k] < 0.0
      ideal[k] += 1.0
    end
    while ideal[k] > 1.0
      ideal[k] -= 1.0
    end
    #println("k: ",k,"  move: ",move,"  new ideal: ",ideal[k])
  end
end

function means( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function attr_vars( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  num_attributes = length(variant_table[1].attributes)
  #println("attr_vars: num_attributes: ",num_attributes)
  ave_vars = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_vars = [ var([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    #println(s," att_vars: ",att_vars)
    ave_vars[i] = mean(att_vars)
    i += 1
  end
  #println("ave_vars: ",ave_vars)
  return ave_vars
end

function count_pops_below_fit_cutoff( mmeans::Vector{Float64}, opt_loss_cutoff::Float64 )
  count = 0
  #for fm = mmeans  # modified 4/4/17
    fm = maximum(mmeans)   # added 4/4/17
    if fm < opt_loss_cutoff
      count += 1
    end
  #end
  count
end

#=
run_evolve(N,num_attributes,num_subpops, ngens,mutation_stddev,num_emmigrants,move_range,move_time_interval,opt_loss_cutoff,uniform_start,min_fit=min_fit,
      linear_fitness=linear_fitness)
tr = temporal_result( T, N, num_attributes, num_subpops, ngens, mutation_stddev, num_emmigrants, move_range, move_time_interval,
          opt_loss_cutoff, horiz_select, min_fit, uniform_start=uniform_start, linear_fitness=linear_fitness, burn_in=burn_in )
function ev_init()
  #include("types.jl")
  include("propsel.jl")
  include("fitness.jl")
  global T = 2
  global N=64
  global num_subpops = 8
  global num_attributes = 4
  global ngens = 200
  global num_emmigrants = 2
  global mutation_stddev = 0.04
  global move_range = 0.1
  global move_time_interval = 5
  global opt_loss_cutoff = 0.31
  global horiz_select = false
  global uniform_start = false
  global min_fit = 0.3
  global linear_fitness=true
  global burn_in = 1.0
end  
=#
