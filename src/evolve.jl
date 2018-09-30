export temporal_result, repeat_evolve, evolve, mutate_meta_pop!, fit_means, fit_means_vars_stddevs, init_meta_pop,
  count_individuals_at_minFit, count_subpops_at_minFit
    
#include("types.jl")
#include("propsel.jl")
#include("fitness.jl") 
#include("mutate_debug.jl")

@doc """ function repeat_evolve( )
Runs evolve()  tr[:num_subpops] times, and averages the results.
See types.jl for the definition of param_type and result_type, and see types.jl for legal keys for this type.
"""
function repeat_evolve( tp::param_type, tr::result_type )  # tp is parameter dictionary, tr is results dictionary
  println("rep_evolve: num_subpops: ",tp[:num_subpops],"  ngens: ",tp[:ngens],"  num_emigrants: ",tp[:num_emigrants],"  horiz_sel: ",tp[:horiz_select],"  mutStddev: ",tp[:mutStddev],
        "  topology: ",tp[:topology],"  num_attributes: ",tp[:num_attributes])
  if tp[:num_trials] == 1
    return evolve( tp, tr )
  end
  int_burn_in = Int(round(tp[:burn_in]*tp[:N]))
  #println("int_burn_in: ",int_burn_in)
  tr_list = result_type[]
  sum_fit_mean = 0.0
  sum_fit_vars = 0.0
  sum_fit_stddev = 0.0
  sum_fit_coef_var = 0.0
  sum_attr_mean = 0.0
  sum_attr_vars = 0.0
  sum_attr_stddev = 0.0
  sum_attr_coef_var = 0.0
  sum_mean_fraction_subpops_below_minFit = 0
  sum_fract_gens_with_all_subpops_below_minFit = 0
  sum_innovations_per_gen_trial = 0
  sum_half_innovations_per_gen_trial = 0
  sum_deleterious_mutations_per_gen_trial = 0
  sum_half_deleterious_mutations_per_gen_trial = 0

  for t = 1:tp[:num_trials]
    tr = evolve( tp, tr ) 
    Base.push!( tr_list, deepcopy(tr) )
    #println("t: ",t," tp[:fitness_variance]: ",tr[:fitness_variance])
    sum_fit_mean += tr[:fitness_mean]
    sum_fit_vars += tr[:fitness_variance]
    sum_fit_stddev += tr[:fitness_stddev]
    sum_fit_coef_var += tr[:fitness_coef_var]
    sum_attr_mean += tr[:attribute_mean]
    sum_attr_vars += tr[:attribute_variance]
    sum_attr_stddev += tr[:attribute_stddev]
    sum_attr_coef_var += tr[:attribute_coef_var]
    sum_innovations_per_gen_trial += tr[:innovations_per_gen_trial]
    sum_half_innovations_per_gen_trial += tr[:half_innovations_per_gen_trial]
    sum_deleterious_mutations_per_gen_trial += tr[:deleterious_mutations_per_gen_trial]
    sum_half_deleterious_mutations_per_gen_trial += tr[:half_deleterious_mutations_per_gen_trial]
    sum_mean_fraction_subpops_below_minFit+= tr[:mean_fraction_subpops_below_minFit]
    sum_fract_gens_with_all_subpops_below_minFit += tr[:fraction_gens_with_all_subpops_below_minFit]
  end
  #println("  tp[:num_trials]: ",tp[:num_trials])
  tr[:mean_fraction_subpops_below_minFit] = sum_mean_fraction_subpops_below_minFit/tp[:num_trials]
  tr[:fraction_gens_with_all_subpops_below_minFit] = sum_fract_gens_with_all_subpops_below_minFit/tp[:num_trials]
  tr[:fitness_mean] = sum_fit_mean/tp[:num_trials]
  tr[:fitness_variance] = sum_fit_vars/tp[:num_trials]
  tr[:fitness_stddev] = sum_fit_stddev/tp[:num_trials]
  tr[:fitness_coef_var] = sum_fit_coef_var/tp[:num_trials]
  tr[:attribute_mean] = sum_attr_vars/tp[:num_trials]
  tr[:attribute_variance] = sum_attr_vars/tp[:num_trials]
  tr[:attribute_stddev] = sum_attr_stddev/tp[:num_trials]
  tr[:attribute_coef_var] = sum_attr_coef_var/tp[:num_trials]
  tr[:innovations_per_gen_trial] = sum_innovations_per_gen_trial/tp[:num_trials]
  tr[:half_innovations_per_gen_trial] = sum_half_innovations_per_gen_trial/tp[:num_trials]
  tr[:deleterious_mutations_per_gen_trial] = sum_deleterious_mutations_per_gen_trial/tp[:num_trials]
  tr[:half_deleterious_mutations_per_gen_trial] = sum_half_deleterious_mutations_per_gen_trial/tp[:num_trials]
  #println( "innov_counts.pos:",tr[:innovations_per_gen_trial], "  innov_counts.half_pos:",tr[:half_innovations_per_gen_trial], 
  #      "  innov_counts.neg:",tr[:deleterious_mutations_per_gen_trial], "  innov_counts.half_neg:",tr[:half_deleterious_mutations_per_gen_trial])
  return tr
end

@doc """ function evolve( )
The main simulation function for temporal.
See types.jl for the definition of param_type and result_type, and see run.jl for legal keys for this type.
"""
function evolve( tp::param_type, tr::result_type ) # tp is parameter dictionary, tr is results dictionary
  #println("function evolve") ###
  #println("num_subpops: ",tp[:num_subpops],"  num_emigrants: ",tp[:num_emigrants],"  horiz_sel: ",tp[:horiz_select],"  mutStddev: ",tp[:mutStddev],"  topology: ",tp[:topology])
  int_burn_in = Int(round(tp[:burn_in]*tp[:N]))
  #println("int_burn_in: ",int_burn_in)
  id = [0]
  # Initialize variables that accumulate over generations
  #att_vars = zeros(tp[:num_subpops])
  cumm_fit_means = zeros(tp[:num_subpops])
  cumm_fit_vars = zeros(tp[:num_subpops])
  cumm_fit_stddevs = zeros(tp[:num_subpops])
  cumm_fit_coef_vars = zeros(tp[:num_subpops])
  cumm_attr_means = zeros(tp[:num_subpops])
  cumm_attr_vars = zeros(tp[:num_subpops])
  cumm_attr_stddevs = zeros(tp[:num_subpops])
  cumm_attr_coef_vars = zeros(tp[:num_subpops])
  cumm_count_subpops_below_minFit = 0
  cumm_count_gens_with_all_subpops_below_minFit = 0
  cumm_count_avg_below_cutoff = 0
  trial_innov_counts = trial_innovation_counts(0.0,0.0,0.0,0.0)
  ideal = fill( tp[:ideal_init], tp[:num_attributes] )
  subpop_size = Int(floor(tp[:N]/tp[:num_subpops]))
  if subpop_size*tp[:num_subpops] != tp[:N]
    println("N:",tp[:N],"  subpop_size: ",subpop_size,"  tp[:num_subpops]: ",tp[:num_subpops],"  prod: ",subpop_size*tp[:num_subpops])
    error("N must equal subpop_size*tp[:num_subpops]")
  end
  vt = Dict{Int64,temporal_variant_type}()
  meta_pop = init_meta_pop( tp, vt, ideal, id )
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  # Iterate over generations
  for g = 1:(tp[:ngens]+int_burn_in)
    if g > int_burn_in && tp[:move_time_interval] > 0 && g % tp[:move_time_interval] == 0
      move_optima( ideal, tp[:move_range] )
      #println("optimum moved:  ideal: ",ideal)
    end
    # mutate and compute gen_innov_counts
    gen_innov_counts = mutate_meta_pop!( meta_pop, vt, ideal, id, tp )  # will also re-evaluate fitness
    #println("after mutate: meta_pop: ",meta_pop)
    #println("gen_innov_counts: ",gen_innov_counts)
    for  j = 1:tp[:num_subpops]
      propsel!( meta_pop[j], vt )  # comment out for fitness test
    end
    #=
    println("after propsel: meta_pop: ",meta_pop)
    for  j = 1:tp[:num_subpops]
      println("j: ",j," ",[ (vt[v].attributes[1], vt[v].fitness) for v in meta_pop[j] ])
    end 
    =#
    mmeans, vvars, stddevs, cfvars = fit_means_vars_stddevs_cfvars( meta_pop, vt )  # just to compute mmeans for horiz_transfer
    horiz_transfer( meta_pop, tp, vt, ideal, mmeans, id, g )
    #println("vt: ",vt)
    if g > int_burn_in  # data collection
      mmeans, vvars, stddevs, cfvars = fit_means_vars_stddevs_cfvars( meta_pop, vt )
      cumm_fit_means += mmeans
      cumm_fit_vars += vvars
      cumm_fit_stddevs += stddevs
      cumm_fit_coef_vars += cfvars
      mmeans, vvars, stddevs, cfvars = attr_means_vars_stddevs_cfvars( meta_pop, vt )
      cumm_attr_means += mmeans
      cumm_attr_vars += vvars
      cumm_attr_stddevs += stddevs
      cumm_attr_coef_vars += cfvars
      # gen_innov_counts was updated by mutate_meta_pop! above.
      trial_innov_counts.pos += gen_innov_counts.pos
      trial_innov_counts.half_pos += gen_innov_counts.half_pos
      trial_innov_counts.neg += gen_innov_counts.neg
      trial_innov_counts.half_neg += gen_innov_counts.half_neg
      count_subpops_below_minFit = count_subpops_at_minFit( meta_pop, vt, tp[:minFit] )
      count_gens_with_all_subpops_below_minFit = (count_subpops_below_minFit == tp[:num_subpops]) ? 1 : 0
      cumm_count_subpops_below_minFit += count_subpops_below_minFit
      cumm_count_gens_with_all_subpops_below_minFit += count_gens_with_all_subpops_below_minFit
      #=  uncomment for fitness test
      for  j = 1:tp[:num_subpops]
        #println("j: ",j," ",[ (vt[v].attributes[1], vt[v].fitness) for v in meta_pop[j] ])
        println("j: ",j," ",[ vt[v].attributes[1] for v in meta_pop[j] ])
        println("variance: ",var([ vt[v].attributes[1] for v in meta_pop[j]]),"  std: ",std([ vt[v].attributes[1] for v in meta_pop[j]]))
      end
      =#
      #println("g: ",g,"  metapop: ",meta_pop)
      #print("g: ",g,"  mmeans: ",mmeans,"  ")
      #println("gen_innov_counts.pos:",gen_innov_counts.pos,"  gen_innov_counts.half_pos:",gen_innov_counts.half_pos,"  gen_innov_counts.neg:",gen_innov_counts.neg,"  gen_innov_counts.half_neg:",gen_innov_counts.half_neg)
      #println( "  tr[:half_deleterious_mutations_per_gen_trial]: ",tr[:half_deleterious_mutations_per_gen_trial])
      #println("   fstdev: ",sqrt(vvars))
      #att_vars = attr_vars( meta_pop, vt )
      #cumm_attr_vars += att_vars
      #println("   astdev: ",sqrt(att_vars))
      #count_gens_with_all_subpops_below_minFit += count_subpops_below_minFit 
      ##print("  count_subpops_below_minFit: ",count_subpops_below_minFit)
      ##println("  count_gens_with_all_subpops_below_minFit: ",count_gens_with_all_subpops_below_minFit)
      #println("cumm_count_subpops_below_minFit: ",cumm_count_subpops_below_minFit)
      #println("cumm_count_gens_with_subpop_below_minFit: ",cumm_count_gens_with_all_subpops_below_minFit)
    end
  end
  # in the following 8 lines, means are taken over subpops
  tr[:fitness_mean] = mean(cumm_fit_means/tp[:ngens])
  tr[:fitness_variance] = mean(cumm_fit_vars/tp[:ngens])
  tr[:fitness_stddev] = mean(cumm_fit_stddevs/tp[:ngens])
  tr[:fitness_coef_var] = mean(cumm_fit_coef_vars/tp[:ngens])
  tr[:attribute_mean] = mean(cumm_attr_means/tp[:ngens])
  tr[:attribute_variance] = mean(cumm_attr_vars/tp[:ngens])
  tr[:attribute_stddev] = mean(cumm_attr_stddevs/tp[:ngens])
  tr[:attribute_coef_var] = mean(cumm_attr_coef_vars/tp[:ngens])
  #println("cumm_fit_means/tp[:ngens]: ",cumm_fit_means/tp[:ngens],"  mean: ",tr[:fitness_mean])
  #println("fit mean: ", tr[:fitness_mean])
  #println("cumm_count_subpops_below_minFit: ",cumm_count_subpops_below_minFit)
  #println("cumm_count_gens_with_subpop_below_minFit: ",cumm_count_gens_with_all_subpops_below_minFit)
  tr[:innovations_per_gen_trial] = trial_innov_counts.pos/tp[:ngens]
  tr[:half_innovations_per_gen_trial] = trial_innov_counts.half_pos/tp[:ngens]
  tr[:deleterious_mutations_per_gen_trial] = trial_innov_counts.neg/tp[:ngens]
  tr[:half_deleterious_mutations_per_gen_trial] = trial_innov_counts.half_neg/tp[:ngens]
  tr[:mean_fraction_subpops_below_minFit] = cumm_count_subpops_below_minFit/tp[:ngens]/tp[:num_subpops]
  tr[:fraction_gens_with_all_subpops_below_minFit] = cumm_count_gens_with_all_subpops_below_minFit/tp[:ngens]
  #println("  innov_counts.pos:",tr[:innovations_per_gen_trial], "  innov_counts.half_pos:",tr[:half_innovations_per_gen_trial], 
  #      "  innov_counts.neg:",tr[:deleterious_mutations_per_gen_trial], "  innov_counts.half_neg:",tr[:half_deleterious_mutations_per_gen_trial])
  #println("fitness variance: ", tr[:fitness_variance],"  attribute variance: ",  tr[:attribute_variance])
  return tr
end

function init_meta_pop( tp::param_type, vt::Dict{Int64,temporal_variant_type}, ideal::Vector{Float64}, id::Vector{Int64} )
  #println("init_meta_pop  meta_pop: ",meta_pop)
  subpop_size = Int(floor(tp[:N]/tp[:num_subpops]))
  if tp[:uniform_start]
    meta_pop = [  fill(1,subpop_size)  for j = 1:tp[:num_subpops] ]
    attr = ideal
    vt[1] = temporal_variant_type( fitness( attr, ideal, minFit=tp[:minFit], linfit_slope=tp[:linfit_slope] ), attr, 0 )
    id[1] += 1
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tp[:num_subpops] ]
    for i = 1:tp[:N]
      attr = rand(tp[:num_attributes])
      vt[i] = temporal_variant_type( fitness( attr, ideal, minFit=tp[:minFit], linfit_slope=tp[:linfit_slope] ), attr, 0 )
    end
    id[1] += tp[:N]
  end
  #println("end init_meta_pop  meta_pop:",meta_pop)
  #println("end init_meta_pop")
  return meta_pop
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

@doc """ count_individuals_at_minFit()
  Counts the number of individuals in a subpop whose fitness is less than or equal to minFit
"""
function count_individuals_at_minFit( subpop::Population, variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
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

@doc """ function count_subpops_at_minFit( )
  Counts the number of subpops where the fitness of all individuals is less than or equal to minFit
"""
function count_subpops_at_minFit( meta_pop::PopList, variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  count = 0
  for j = 1:num_subpops
    if count_individuals_at_minFit( meta_pop[j], variant_table, minFit ) == subpop_size
      count += 1
    end
  end
  ##println("count_subpops_at_minFit: ",count,"  subpop_size: ",subpop_size,"  minFit: ",minFit)
  return count
end


#=
run_evolve(N,num_attributes,num_subpops, ngens,mutStddev,num_emigrants,move_range,move_time_interval,minFit,uniform_start,minFit=minFit)
tr = temporal_result( T, N, num_attributes, num_subpops, ngens, mutStddev, num_emigrants, move_range, move_time_interval,
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
  global num_emigrants = 2
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
