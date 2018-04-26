export repeat_evolve_until_dead, evolve_until_dead

@doc """ function repeat_evolve_until_dead()
  Run repeated trails of evolve_until_dead()
"""
function repeat_evolve_until_dead( paramd::param_type, resultd::result_type )
  println("repeat_evolve_until_dead: N: ",paramd[:N],"  num_attributes: ",paramd[:num_attributes]," num_subpops: ",paramd[:num_subpops],"  num_emigrants: ",paramd[:num_emigrants])
  if paramd[:num_trials] == 1
    return evolve_until_dead( paramd, resultd )
  end
  resultd_list = result_type[]
  sum_generational_lifetime = 0.0
  sum_move_update_lifetime = 0.0
  sum_gen_limit_count = 0
  for t = 1:paramd[:num_trials]
    resultd = evolve_until_dead( paramd, deepcopy(resultd) )
    Base.push!( resultd_list, deepcopy(resultd) )   # TODO:  what is this line doing?
    sum_generational_lifetime += resultd[:generational_lifetime] 
    sum_move_update_lifetime += resultd[:move_update_lifetime] 
    sum_gen_limit_count += resultd[:gen_limit_reached_count]
  end
  resultd[:generational_lifetime] = sum_generational_lifetime/paramd[:num_trials]
  resultd[:move_update_lifetime] = sum_move_update_lifetime/paramd[:num_trials]
  resultd[:gen_limit_reached_count] = sum_gen_limit_count
  return resultd
end

@doc """ function evolve_until_dead( )
A subpop is "dead" if all members have fitness minFit.
Evolve generations until all subpops are dead, either on a single generation update, or on a move optimum update
Doesn't keep track of the other statistics.
See types.jl for the definition of param_type and result_type, and for the definition of the fields of this type.
"""
function evolve_until_dead( paramd::param_type, resultd::result_type )
  #println("function evolve_until_dead")  ###
  #println("num_subpops: ",paramd[:num_subpops],"  num_emmigrants: ",paramd[:num_emigrants],"  horiz_sel: ",paramd[:horiz_select],"  mutStddev: ",paramd[:mutStddev],"  topology: ",paramd[:topology])
  #println("burn_in: ",paramd[:burn_in])
  int_burn_in = Int(round(paramd[:burn_in]*paramd[:N]))
  id = [0]
  ideal = fill( paramd[:ideal_init], paramd[:num_attributes] )
  subpop_size = Int(floor(paramd[:N]/paramd[:num_subpops]))
  if subpop_size*paramd[:num_subpops] != paramd[:N]
    error("N must equal subpop_size*paramd[:num_subpops]")
  end
  vt = Dict{Int64,variant_type}()
  meta_pop = init_meta_pop( paramd, vt, ideal, id )
  #subpop_alive = fill(true,paramd[:num_subpops])  # Should be used with uniform start
  #prev_subpop_alive = fill(true,paramd[:num_subpops])  # Should be used with uniform start
  sbp = subpop_properties_init( paramd[:num_subpops])
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  resultd[:generational_lifetime]=0
  resultd[:move_update_lifetime]=0
  resultd[:gen_limit_reached_count]=0
  g = 1
  #while (any(sbp.generational_subpop_alive) || any(sbp.prev_subpop_alive)) && g < paramd[:ngens]+int_burn_in
  while (resultd[:generational_lifetime]==0 || resultd[:move_update_lifetime]==0 ) && g < paramd[:ngens]+int_burn_in
    if paramd[:move_time_interval] > 0 && g > int_burn_in && g % paramd[:move_time_interval] == 0
      #println(" g: ",g,"  #optimum to be moved   means: ", fmeans( meta_pop, vt ))
      #println("B g: ",g,"  sbp.prev_subpop_alive: ",sbp.prev_subpop_alive )
      #println("B g: ",g,"  sbp.current_subpop_alive: ",sbp.current_subpop_alive )
      move_optima( ideal, paramd[:move_range] )
      #println(" g: ",g,"  #optimum just  moved   means: ", fmeans( meta_pop, vt ))
      subpop_alive_opt_move_update( sbp, meta_pop,  vt, paramd[:minFit] )
      #println("A g: ",g,"  sbp.prev_subpop_alive: ",sbp.prev_subpop_alive )
      #println("A g: ",g,"  sbp.current_subpop_alive: ",sbp.current_subpop_alive )
    end
    mutate_meta_pop!( meta_pop, vt, ideal, id, paramd )  # will also re-evaluate fitness
    mmeans = fmeans( meta_pop, vt )
    #println("  g: ",g," after mutate means: ",mmeans)
    for  j = 1:paramd[:num_subpops]
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )  # comment out for fitness test
    end
    mmeans = fmeans( meta_pop, vt )
    #println("  g: ",g," before horiz means: ",mmeans)
    horiz_transfer( meta_pop, paramd, vt, ideal, mmeans, id, g )
    mmeans, vvars = means_vars( meta_pop, vt )
    #println("  g: ",g," after horiz: ",mmeans)
    #println("vt: ",vt)
    subpop_alive_gen_update( sbp, meta_pop, vt, paramd[:minFit] )
    #println("  g: ",g,"  sbp.generational_subpop_alive: ",sbp.generational_subpop_alive )
    #println("  g: ",g,"  sbp.prev_subpop_alive: ",sbp.prev_subpop_alive )
    #println("  g: ",g,"  sbp.current_subpop_alive: ",sbp.current_subpop_alive )
    if resultd[:generational_lifetime]==0 && !any(sbp.generational_subpop_alive)
      #println("g: ",g," All subpops died (single generation update).")
      resultd[:generational_lifetime] = g
    elseif resultd[:move_update_lifetime]==0 && !any(sbp.prev_subpop_alive)
      #println("g: ",g," All subpops died (move optimum update).")
      resultd[:move_update_lifetime] = g
    end
    g += 1
  end
  if g == paramd[:ngens]+int_burn_in
    #println("evolve_until_dead ran for ",g," generations when the generation limit was reached")
    resultd[:gen_limit_reached_count] += 1
    if resultd[:generational_lifetime] > 0
      #println("all subpops died at generation ",resultd[:generational_lifetime],"  (generational update).")
    else
      resultd[:generational_lifetime] = g
    end
    if resultd[:move_update_lifetime] > 0
      #println("all subpops died at generation ",resultd[:move_update_lifetime],"(move optimum update).")
    else
      resultd[:move_update_lifetime] = g
    end
  #else 
    #println("all subpops died at generation ",resultd[:generational_lifetime],"  (genupdate).")
    #println("all subpops died at generation ",resultd[:move_update_lifetime],"  (move optimum update).")
  end
  return resultd
end

@doc """ function subpop_properties_init()
   Initializes 3 lists of Bools that describe the 'alive' properties of the subpops
"""
function subpop_properties_init( num_subpops::Int64 )
  generational_subpop_alive = fill(true,num_subpops)
  prev_subpop_alive = fill(true,num_subpops)
  current_subpop_alive = fill(false,num_subpops)
  return subpop_properties(generational_subpop_alive,prev_subpop_alive,current_subpop_alive)
end

function subpop_alive_gen_update( sbp::subpop_properties, meta_pop::PopList,  variant_table::Dict{Int64,variant_type}, minFit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  for j = 1:num_subpops
    count_indivs_below_minfit = count_individuals_below_minfit( meta_pop[j], variant_table, minFit )
    if sbp.generational_subpop_alive[j] && count_indivs_below_minfit == subpop_size
      sbp.generational_subpop_alive[j] = false  # only set to false if previously true and all individuals have fitness <= minFit
    end
    if !sbp.current_subpop_alive[j] && count_indivs_below_minfit < subpop_size
      sbp.current_subpop_alive[j] = true  # set to true if previously false and some individual has fitness greater than minFit
    end
  end
end

function subpop_alive_opt_move_update( sbp::subpop_properties, meta_pop::PopList,  variant_table::Dict{Int64,variant_type}, minFit::Float64 )
  if !any(sbp.prev_subpop_alive) 
    return
  end
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  for j = 1:num_subpops
    if sbp.prev_subpop_alive[j] && !sbp.current_subpop_alive[j]
      sbp.prev_subpop_alive[j] = false  # only set to false if previously true and all individuals have fitness less than minFit
    end
    sbp.current_subpop_alive[j] = false
  end
end

@doc """ function update_subpop_alive!()
  A subpop is alive if it satisfies two conditions: it was alive in the previous generation, and it has at least one individual whose
     fitness is strictly greater than minFit.
"""
function update_subpop_alive!( subpop_alive::Vector{Bool}, prev_subpop_alive::Vector{Bool}, meta_pop::PopList, variant_table::Dict{Int64,variant_type}, minFit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  for j = 1:num_subpops
    if subpop_alive[j] && count_individuals_below_minfit( meta_pop[j], variant_table, minFit ) == subpop_size
      subpop_alive[j] = false  # only set to false if previously true and all individuals have fitness less than minFit
    end
  end
end

