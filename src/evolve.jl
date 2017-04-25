export temporal_result, repeat_evolve, evolve, mutate_meta_pop!, means,
  count_individuals_below_cutoff, count_subpops_below_minfit
    
#include("types.jl")
#include("propsel.jl")
#include("fitness.jl") 

function temporal_result( num_trials::Int64, N::Int64, num_attributes::Int64, num_subpops::Int64, ngens::Int64, mutation_stddev::Float64, num_emmigrants::Int64,
    move_range::Float64, move_time_interval::Int64, horiz_select::Bool=false, min_fit::Float64=0.0; topology::String="circular", 
    uniform_start::Bool=false, linear_fitness::Bool=false, linfit_slope::Float64=1.0, burn_in::Float64=1.0 )
  ideal_init = 0.5
  return temporal_result_type( num_trials, N, num_subpops, num_emmigrants, num_attributes, ngens, burn_in, uniform_start, horiz_select, mutation_stddev,
      ideal_init, move_range, move_time_interval, min_fit, linear_fitness, linfit_slope, topology, 0.0, 0.0, 0.0, 0.0, 0.0 )
end

@doc """ function repeat_evolve( )
Runs evolve()  tr.num_subpops times, and averages the results.
"""
function repeat_evolve( tr::temporal_result_type )
  println("rep_evolve: num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutation stddev: ",tr.mutation_stddev,
        "  topology: ",tr.topology,"  linfit_slope: ",tr.linfit_slope)
  if tr.num_trials == 1
    return evolve( tr )
  end
  tr_list = temporal_result_type[]
  sum_mean = 0.0
  sum_vars = 0.0
  sum_attr_vars = 0.0
  sum_mean_fraction_subpops_below_min_fit = 0
  sum_fract_gens_with_all_subpops_below_min_fit = 0

  for t = 1:tr.num_trials
    tr = evolve( tr ) 
    Base.push!( tr_list, deepcopy(tr) )
    sum_mean_fraction_subpops_below_min_fit+= tr.mean_fraction_subpops_below_min_fit
    sum_fract_gens_with_all_subpops_below_min_fit += tr.fraction_gens_with_all_subpops_below_min_fit
    #println("t: ",t," tr.fitness_variance: ",tr.fitness_variance)
    sum_mean += tr.fitness_mean
    sum_vars += tr.fitness_variance
    sum_attr_vars += tr.attribute_variance
  end
  tr.mean_fraction_subpops_below_min_fit = sum_mean_fraction_subpops_below_min_fit/tr.num_trials
  tr.fraction_gens_with_all_subpops_below_min_fit = sum_fract_gens_with_all_subpops_below_min_fit/tr.num_trials
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
  att_vars = zeros(tr.num_subpops)
  cumm_means = zeros(tr.num_subpops)
  cumm_vars = zeros(tr.num_subpops)
  cumm_attr_vars = zeros(tr.num_subpops)
  cumm_count_subpops_below_min_fit = 0
  cumm_count_gens_with_all_subpops_below_min_fit = 0
  cumm_count_avg_below_cutoff = 0
  ideal = fill( tr.ideal_init, tr.num_attributes )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  if subpop_size*tr.num_subpops != tr.N
    error("N must equal subpop_size*tr.num_subpops")
  end
  vt = Dict{Int64,variant_type}()
  meta_pop = init_meta_pop( tr, vt, ideal, id )
  subpop_live = fill(true,tr.num_subpops)
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  for g = 1:(tr.ngens+int_burn_in)
    if g > int_burn_in && g % tr.move_time_interval == 0
      move_optima( ideal, tr.move_range )
      #println("optimum moved")
      #mmeans, vvars = means_vars( meta_pop, vt )
      #println("g: ",g,"  metapop: ",meta_pop)
      #println("g: ",g,"  mmeans: ",mmeans,"  ")
    end
    mutate_meta_pop!( meta_pop, vt, ideal, id, tr )  # will also re-evaluate fitness
    
    for  j = 1:tr.num_subpops
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )  # comment out for fitness test
    end
    
    mmeans = fmeans( meta_pop, vt )
    horiz_transfer( meta_pop, tr, vt, ideal, mmeans, id, g )
    mmeans, vvars = means_vars( meta_pop, vt )
    #println("g: ",g,"  mmeans: ",mmeans)
    #println("vvars: ",vvars)
    #println("vt: ",vt)
    if g > int_burn_in  # data collection
      #(mmeans, vvars) = means( meta_pop, vt )
      cumm_means += mmeans
      cumm_vars += vvars
      #=  uncomment for fitness test
      for  j = 1:tr.num_subpops
        #println("j: ",j," ",[ (vt[v].attributes[1], vt[v].fitness) for v in meta_pop[j] ])
        println("j: ",j," ",[ vt[v].attributes[1] for v in meta_pop[j] ])
        println("variance: ",var([ vt[v].attributes[1] for v in meta_pop[j]]),"  std: ",std([ vt[v].attributes[1] for v in meta_pop[j]]))
      end
      =#
      ##println("g: ",g,"  metapop: ",meta_pop)
      ##println("g: ",g,"  mmeans: ",mmeans,"  ")
      #print("   fstdev: ",sqrt(vvars))
      att_vars = attr_vars( meta_pop, vt )
      cumm_attr_vars += att_vars
      #println("   astdev: ",sqrt(att_vars))
      count_subpops_below_min_fit = count_subpops_below_minfit( meta_pop, vt, tr.min_fit )
      count_gens_with_all_subpops_below_min_fit = (count_subpops_below_min_fit == tr.num_subpops) ? 1 : 0
      #count_gens_with_all_subpops_below_min_fit += count_subpops_below_min_fit 
      ##print("  count_subpops_below_min_fit: ",count_subpops_below_min_fit)
      ##println("  count_gens_with_all_subpops_below_min_fit: ",count_gens_with_all_subpops_below_min_fit)
      cumm_count_subpops_below_min_fit += count_subpops_below_min_fit
      cumm_count_gens_with_all_subpops_below_min_fit += count_gens_with_all_subpops_below_min_fit
      #println("cumm_count_subpops_below_min_fit: ",cumm_count_subpops_below_min_fit)
      #println("cumm_count_gens_with_subpop_below_min_fit: ",cumm_count_gens_with_all_subpops_below_min_fit)
    end
  end
  tr.mean_fraction_subpops_below_min_fit = cumm_count_subpops_below_min_fit/tr.ngens/tr.num_subpops
  tr.fraction_gens_with_all_subpops_below_min_fit = cumm_count_gens_with_all_subpops_below_min_fit/tr.ngens  
  tr.fitness_mean = mean(cumm_means/tr.ngens)
  #println("cumm_means/tr.ngens: ",cumm_means/tr.ngens,"  mean: ",tr.fitness_mean)
  #println("cumm_count_subpops_below_min_fit: ",cumm_count_subpops_below_min_fit)
  #println("cumm_count_gens_with_subpop_below_min_fit: ",cumm_count_gens_with_all_subpops_below_min_fit)
  tr.fitness_variance = mean(cumm_vars/tr.ngens)
  tr.attribute_variance = mean(cumm_attr_vars/tr.ngens)
  return tr
end

function init_meta_pop( tr::temporal_result_type, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64} )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  if tr.uniform_start
    meta_pop = [  fill(1,subpop_size)  for j = 1:tr.num_subpops ]
    attr = ideal
    vt[1] = variant_type( fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness, linfit_slope=tr.linfit_slope ), attr )
    id[1] += 1
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tr.num_subpops ]
    for i = 1:tr.N
      attr = rand(tr.num_attributes)
      vt[i] = variant_type( fitness( attr, ideal, min_fit=tr.min_fit, linear_fitness=tr.linear_fitness, linfit_slope=tr.linfit_slope ), attr )
    end
    id[1] += tr.N
  end
  return meta_pop
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
      vt[meta_pop[j][i]].fitness = fitness( vt[meta_pop[j][i]].attributes, ideal, min_fit=tr.min_fit, 
          linear_fitness=tr.linear_fitness, linfit_slope=tr.linfit_slope )
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

function means_vars( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function fmeans( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  return means
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

@doc """ count_individuals_below_minfit()
  Counts the number of individuals in a subpop whose fitness is less than or equal to min_fit
"""
function count_individuals_below_minfit( subpop::Population, variant_table::Dict{Int64,variant_type}, min_fit::Float64 )
  count = 0
  for v in subpop
    if variant_table[v].fitness <= min_fit + 10.0*eps()
      count += 1
    end
    #print(v,":",variant_table[v].fitness,"; ")
  end
  ##println("   count indivs below minfit: ",count)
  return count
end

@doc """ function count_subpops_below_minfit( )
  Counts the number of subpops where the fitness of all individuals is less than or equal to min_fit
"""
function count_subpops_below_minfit( meta_pop::PopList, variant_table::Dict{Int64,variant_type}, min_fit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  count = 0
  for j = 1:num_subpops
    if count_individuals_below_minfit( meta_pop[j], variant_table, min_fit ) == subpop_size
      count += 1
    end
  end
  ##println("count_subpops_below_minfit: ",count,"  subpop_size: ",subpop_size,"  min_fit: ",min_fit)
  return count
end

#=
run_evolve(N,num_attributes,num_subpops, ngens,mutation_stddev,num_emmigrants,move_range,move_time_interval,min_fit,uniform_start,min_fit=min_fit,
      linear_fitness=linear_fitness)
tr = temporal_result( T, N, num_attributes, num_subpops, ngens, mutation_stddev, num_emmigrants, move_range, move_time_interval,
          min_fit, horiz_select, min_fit, uniform_start=uniform_start, linear_fitness=linear_fitness, burn_in=burn_in )
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
  global min_fit = 0.31
  global horiz_select = false
  global uniform_start = false
  global min_fit = 0.3
  global linear_fitness=true
  global burn_in = 1.0
end  
=#
