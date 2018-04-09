export temporal_result, repeat_evolve, evolve, mutate_meta_pop!, fmeans, means_vars, init_meta_pop,
  count_individuals_below_minfit, count_subpops_below_minfit
    
#include("types.jl")
#include("propsel.jl")
#include("fitness.jl") 
#include("mutate_debug.jl")

@doc """ function repeat_evolve( )
Runs evolve()  tr.num_subpops times, and averages the results.
"""
function repeat_evolve( tr::temporal_result_type )
  #println("rep_evolve: num_subpops: ",tr.num_subpops,"  ngens: ",tr.ngens,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutStddev: ",tr.mutStddev,
  #      "  topology: ",tr.topology,"  num_attributes: ",tr.num_attributes)
  if tr.num_trials == 1
    return evolve( tr )
  end
  int_burn_in = Int(round(tr.burn_in*tr.N))
  #println("int_burn_in: ",int_burn_in)
  tr_list = temporal_result_type[]
  sum_mean = 0.0
  sum_vars = 0.0
  sum_attr_vars = 0.0
  sum_mean_fraction_subpops_below_minFit = 0
  sum_fract_gens_with_all_subpops_below_minFit = 0
  sum_innovations_per_gen_trial = 0
  sum_half_innovations_per_gen_trial = 0
  sum_deleterious_mutations_per_gen_trial = 0
  sum_half_deleterious_mutations_per_gen_trial = 0

  for t = 1:tr.num_trials
    tr = evolve( tr ) 
    Base.push!( tr_list, deepcopy(tr) )
    sum_mean_fraction_subpops_below_minFit+= tr.mean_fraction_subpops_below_minFit
    sum_fract_gens_with_all_subpops_below_minFit += tr.fraction_gens_with_all_subpops_below_minFit
    #println("t: ",t," tr.fitness_variance: ",tr.fitness_variance)
    sum_mean += tr.fitness_mean
    sum_vars += tr.fitness_variance
    sum_attr_vars += tr.attribute_variance
    sum_innovations_per_gen_trial += tr.innovations_per_gen_trial
    sum_half_innovations_per_gen_trial += tr.half_innovations_per_gen_trial
    sum_deleterious_mutations_per_gen_trial += tr.deleterious_mutations_per_gen_trial
    sum_half_deleterious_mutations_per_gen_trial += tr.half_deleterious_mutations_per_gen_trial
  end
  tr.mean_fraction_subpops_below_minFit = sum_mean_fraction_subpops_below_minFit/tr.num_trials
  tr.fraction_gens_with_all_subpops_below_minFit = sum_fract_gens_with_all_subpops_below_minFit/tr.num_trials
  tr.fitness_mean = sum_mean/tr.num_trials
  tr.fitness_variance = sum_vars/tr.num_trials
  #println("  tr.fitness_variance: ",tr.fitness_variance)
  #println("  tr.num_trials: ",tr.num_trials)
  tr.attribute_variance = sum_attr_vars/tr.num_trials
  tr.innovations_per_gen_trial = sum_innovations_per_gen_trial/tr.num_trials
  tr.half_innovations_per_gen_trial = sum_half_innovations_per_gen_trial/tr.num_trials
  tr.deleterious_mutations_per_gen_trial = sum_deleterious_mutations_per_gen_trial/tr.num_trials
  tr.half_deleterious_mutations_per_gen_trial = sum_half_deleterious_mutations_per_gen_trial/tr.num_trials
  #println( "innov_counts.pos:",tr.innovations_per_gen_trial, "  innov_counts.half_pos:",tr.half_innovations_per_gen_trial, 
  #      "  innov_counts.neg:",tr.deleterious_mutations_per_gen_trial, "  innov_counts.half_neg:",tr.half_deleterious_mutations_per_gen_trial)
  #return tr, tr_list
  return tr
end

@doc """ function evolve( )
The main simulation function for temporal.
See types.jl for the definition of temporal_result_type, and for the definition of the fields of this type.
"""
function evolve( tr::temporal_result_type )
  #println("function evolve") ###
  #println("num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutStddev: ",tr.mutStddev,"  topology: ",tr.topology)
  int_burn_in = Int(round(tr.burn_in*tr.N))
  #println("int_burn_in: ",int_burn_in)
  id = [0]
  att_vars = zeros(tr.num_subpops)
  cumm_means = zeros(tr.num_subpops)
  cumm_vars = zeros(tr.num_subpops)
  cumm_attr_vars = zeros(tr.num_subpops)
  cumm_count_subpops_below_minFit = 0
  cumm_count_gens_with_all_subpops_below_minFit = 0
  cumm_count_avg_below_cutoff = 0
  ideal = fill( tr.ideal_init, tr.num_attributes )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  if subpop_size*tr.num_subpops != tr.N
    println("N:",tr.N,"  subpop_size: ",subpop_size,"  tr.num_subpops: ",tr.num_subpops,"  prod: ",subpop_size*tr.num_subpops)
    error("N must equal subpop_size*tr.num_subpops")
  end
  trial_innov_counts = trial_innovation_counts(0.0,0.0,0.0,0.0)
  vt = Dict{Int64,variant_type}()
  meta_pop = init_meta_pop( tr, vt, ideal, id )
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  for g = 1:(tr.ngens+int_burn_in)
    if g > int_burn_in && tr.move_time_interval > 0 && g % tr.move_time_interval == 0
      move_optima( ideal, tr.move_range )
      #println("optimum moved")
      mmeans, vvars = means_vars( meta_pop, vt )
      #println("g: ",g,"  metapop: ",meta_pop)
      #println("g: ",g,"  mmeans: ",mmeans,"  ")
    end
    gen_innov_counts = mutate_meta_pop!( meta_pop, vt, ideal, id, tr )  # will also re-evaluate fitness
    #println("gen_innov_counts: ",gen_innov_counts)
    
    for  j = 1:tr.num_subpops
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )  # comment out for fitness test
    end
    
    mmeans = fmeans( meta_pop, vt )
    #println("g: ",g,"  mmeans: ",mmeans,"  ")
    horiz_transfer( meta_pop, tr, vt, ideal, mmeans, id, g )
    mmeans, vvars = means_vars( meta_pop, vt )
    #println("g: ",g,"  mmeans: ",mmeans,"  ")
    #println("g: ",g,"  mmeans: ",mmeans)
    #println("vvars: ",vvars)
    #println("vt: ",vt)
    if g > int_burn_in  # data collection
      #(mmeans, vvars) = means_vars( meta_pop, vt )
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
      trial_innov_counts.pos += gen_innov_counts.pos
      trial_innov_counts.half_pos += gen_innov_counts.half_pos
      trial_innov_counts.neg += gen_innov_counts.neg
      trial_innov_counts.half_neg += gen_innov_counts.half_neg
      #print("g: ",g,"  mmeans: ",mmeans,"  ")
      #println("gen_innov_counts.pos:",gen_innov_counts.pos,"  gen_innov_counts.half_pos:",gen_innov_counts.half_pos,"  gen_innov_counts.neg:",gen_innov_counts.neg,"  gen_innov_counts.half_neg:",gen_innov_counts.half_neg)
      #println( "  tr.half_deleterious_mutations_per_gen_trial: ",tr.half_deleterious_mutations_per_gen_trial)
      #print("   fstdev: ",sqrt(vvars))
      att_vars = attr_vars( meta_pop, vt )
      cumm_attr_vars += att_vars
      #println("   astdev: ",sqrt(att_vars))
      count_subpops_below_minFit = count_subpops_below_minfit( meta_pop, vt, tr.minFit )
      count_gens_with_all_subpops_below_minFit = (count_subpops_below_minFit == tr.num_subpops) ? 1 : 0
      #count_gens_with_all_subpops_below_minFit += count_subpops_below_minFit 
      ##print("  count_subpops_below_minFit: ",count_subpops_below_minFit)
      ##println("  count_gens_with_all_subpops_below_minFit: ",count_gens_with_all_subpops_below_minFit)
      cumm_count_subpops_below_minFit += count_subpops_below_minFit
      cumm_count_gens_with_all_subpops_below_minFit += count_gens_with_all_subpops_below_minFit
      #println("cumm_count_subpops_below_minFit: ",cumm_count_subpops_below_minFit)
      #println("cumm_count_gens_with_subpop_below_minFit: ",cumm_count_gens_with_all_subpops_below_minFit)
    end
  end
  tr.mean_fraction_subpops_below_minFit = cumm_count_subpops_below_minFit/tr.ngens/tr.num_subpops
  tr.fraction_gens_with_all_subpops_below_minFit = cumm_count_gens_with_all_subpops_below_minFit/tr.ngens  
  tr.fitness_mean = mean(cumm_means/tr.ngens)
  #println("cumm_means/tr.ngens: ",cumm_means/tr.ngens,"  mean: ",tr.fitness_mean)
  #println("cumm_count_subpops_below_minFit: ",cumm_count_subpops_below_minFit)
  #println("cumm_count_gens_with_subpop_below_minFit: ",cumm_count_gens_with_all_subpops_below_minFit)
  tr.innovations_per_gen_trial = trial_innov_counts.pos/tr.ngens
  tr.half_innovations_per_gen_trial = trial_innov_counts.half_pos/tr.ngens
  tr.deleterious_mutations_per_gen_trial = trial_innov_counts.neg/tr.ngens
  tr.half_deleterious_mutations_per_gen_trial = trial_innov_counts.half_neg/tr.ngens
  #println("  innov_counts.pos:",tr.innovations_per_gen_trial, "  innov_counts.half_pos:",tr.half_innovations_per_gen_trial, 
  #      "  innov_counts.neg:",tr.deleterious_mutations_per_gen_trial, "  innov_counts.half_neg:",tr.half_deleterious_mutations_per_gen_trial)
  tr.fitness_variance = mean(cumm_vars/tr.ngens)
  tr.attribute_variance = mean(cumm_attr_vars/tr.ngens)
  return tr
end

function init_meta_pop( tr::temporal_result_type, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64} )
  #println("start init_meta_pop")
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  if tr.uniform_start
    meta_pop = [  fill(1,subpop_size)  for j = 1:tr.num_subpops ]
    attr = ideal
    vt[1] = variant_type( fitness( attr, ideal, minFit=tr.minFit, linfit_slope=tr.linfit_slope ), attr )
    id[1] += 1
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tr.num_subpops ]
    for i = 1:tr.N
      attr = rand(tr.num_attributes)
      vt[i] = variant_type( fitness( attr, ideal, minFit=tr.minFit, linfit_slope=tr.linfit_slope ), attr )
    end
    id[1] += tr.N
  end
  #println("end init_meta_pop  meta_pop:",meta_pop)
  #println("end init_meta_pop")
  return meta_pop
end

function mutate_subpop!( subpop::Population, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, 
    id::Vector{Int64}, tr::temporal_result_type, innov_counts::generational_innovation_counts  )
  subpop_size = length(subpop)
  v_list = [vt[subpop[i]] for i = 1:subpop_size ]
  #innov_counts = generational_innovation_counts(0,0,0,0)
  fitnesses = zeros(Float64,subpop_size)
  for i = 1:subpop_size
    #=
    if id[1] == 480
      println("B i: ",i," subpop[i]: ",subpop[i],"  v_list[i]: ",v_list[i])
    end
    =#
    subpop[i] = id[1]
    id[1] += 1
    vt[subpop[i]] = mutate_variant( v_list[i], tr.mutStddev )
    fitnesses[i] = fitness( vt[subpop[i]].attributes, ideal, minFit=tr.minFit, linfit_slope=tr.linfit_slope )
    vt[subpop[i]].fitness = fitnesses[i]
  end
  mmean = mean(fitnesses)
  num_positive_select = count(x->x>mmean+mmean/subpop_size,fitnesses)
  num_half_positive_select = count(x->x>mmean+0.5*mmean/subpop_size,fitnesses)
  num_negative_select = count(x->x<mmean-1.0*mmean/subpop_size,fitnesses)
  num_half_negative_select = count(x->x<mmean-0.5*mmean/subpop_size,fitnesses)
  #num_half_negative_select = length(fitnesses)
  #println( "num_pos: ",num_positive_select, "  num_half_pos: ",num_half_positive_select, "  num_neg: ",num_negative_select, "  num_half_neg: ",num_half_negative_select)
  #println( "  num_half_neg: ",num_half_negative_select)
  innov_counts.pos += num_positive_select
  innov_counts.half_pos += num_half_positive_select
  innov_counts.neg += num_negative_select
  innov_counts.half_neg += num_half_negative_select
end

@doc """ function mutate_meta_pop!( )
  Mutates each individual of each subpop of meta_pop.  
  Mutation is done by applying the function mutate_variant().
  The mutated individual has a new id, and thus is a different individual.
  For each individual, determines whether the new fitness corresponds to a selection coefficient is in the nearly neutral range.
  The selection coefficient of an individual i of a subpop is defined as  s(i) = fitness(i)/mean_fitness  where mean_fitness is the mean over the subpop.
  Individual i is advantageous if  s(i) > 1 + 1/subpop_size and is disadvantageous if s(i) < 1 - 1/subpopsize.
  We think of advantageous individuals as innovations.
  We also define individual i to be half-advantageous if  s(i) > 1+0.5/subpop_size and half-disadvantagous if s(i) < 1-0.5/subpop_size.
"""
function mutate_meta_pop!( meta_pop::PopList, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64}, tr::temporal_result_type  )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  #println("num_subpops: ",num_subpops,"  subpop_size: ",subpop_size)
  #v_lists = [ [ vt[meta_pop[j][i]] for i = 1:subpop_size ] for j = 1:tr.num_subpops]
  #v_lists = [ [ vt[meta_pop[j][i]] for i = 1:subpop_size ] for j = 1:num_subpops ]
  gen_innov_counts = generational_innovation_counts(0,0,0,0)
  for j = 1:num_subpops
    mutate_subpop!(meta_pop[j],vt,ideal,id,tr,gen_innov_counts)
  end
  return gen_innov_counts
end

function mutate_variant( v::variant_type, mutStddev::Float64 )
  new_v = deepcopy(v)
  new_v.attributes = mutate_attributes( new_v.attributes, mutStddev )
  new_v
end

@doc """ mutate_attributes( )
  Add a normally distributed (with standard deviation mutStddev) random value to each attribute.
  If the result is outside of the unit interval [0,1], wrap around.
"""
function mutate_attributes( attributes::Vector{Float64}, mutStddev::Float64 )
  #println("mutate attributes  mutStddev: ",mutStddev,"  attributes: ",attributes)
  new_attributes = deepcopy(attributes)
  for i = 1:length(attributes)
    #println("B new_attributes[",i,"]: ",new_attributes[i])
    new_attributes[i] += +mutStddev*randn()
    while new_attributes[i] < 0
        new_attributes[i] += 1.0
        #println("wrapped up: ",new_attributes[i])
    end
    while new_attributes[i] > 1.0
        new_attributes[i] -= 1.0
        #println("wrapped down: ",new_attributes[i])
    end
    if new_attributes[i] < 0.0 || new_attributes[i] > 1.0
      error("attribute not wrapped")
    end
    new_attributes[i] = min(1.0,max(0.0,new_attributes[i]))
    #println("A new_attributes[",i,"]: ",new_attributes[i])
  end
  #println("new_attributes: ",new_attributes)
  return new_attributes
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

@doc """ means_vars()
means[s]  is the mean fitness of the individuals of subpopulation s
vars[s]   is the variance of the fitnesses of the individuals of subpopulation s
"""
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

@doc """ attr_vars()
Returns ave_vars, where ave_vars[s] is the mean (over attributes) of the variance of each attribute, where the
  variance is taken over the members of subpopulation s.
"""
function attr_vars( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  num_attributes = length(variant_table[1].attributes)
  #println("attr_vars: num_attributes: ",num_attributes)
  ave_vars = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    # att_vars[i] is variance (taken over the individuals of subpopulation s) of attribute i
    att_vars = [ var([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    #println(s," att_vars: ",att_vars)
    ave_vars[i] = mean(att_vars)   # mean over attributes
    i += 1
  end
  #println("ave_vars: ",ave_vars)
  return ave_vars
end

@doc """ count_individuals_below_minfit()
  Counts the number of individuals in a subpop whose fitness is less than or equal to minFit
"""
function count_individuals_below_minfit( subpop::Population, variant_table::Dict{Int64,variant_type}, minFit::Float64 )
  count = 0
  for v in subpop
    if variant_table[v].fitness <= minFit + 10.0*eps()
      count += 1
    end
    #print(v,":",variant_table[v].fitness,"; ")
  end
  ##println("   count indivs below minfit: ",count)
  return count
end

@doc """ function count_subpops_below_minfit( )
  Counts the number of subpops where the fitness of all individuals is less than or equal to minFit
"""
function count_subpops_below_minfit( meta_pop::PopList, variant_table::Dict{Int64,variant_type}, minFit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  count = 0
  for j = 1:num_subpops
    if count_individuals_below_minfit( meta_pop[j], variant_table, minFit ) == subpop_size
      count += 1
    end
  end
  ##println("count_subpops_below_minfit: ",count,"  subpop_size: ",subpop_size,"  minFit: ",minFit)
  return count
end


#=
run_evolve(N,num_attributes,num_subpops, ngens,mutStddev,num_emmigrants,move_range,move_time_interval,minFit,uniform_start,minFit=minFit)
tr = temporal_result( T, N, num_attributes, num_subpops, ngens, mutStddev, num_emmigrants, move_range, move_time_interval,
          minFit, horiz_select, minFit, uniform_start=uniform_start,  burn_in=burn_in )
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
  global mutStddev = 0.04
  global move_range = 0.1
  global move_time_interval = 5
  global minFit = 0.31
  global horiz_select = false
  global uniform_start = false
  global minFit = 0.3
  global burn_in = 1.0
end  
=#
