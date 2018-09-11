# Convert to new method of computing generational and move-update lifetimes  6/28/18
# This version still contains much debugging code and commented out code for previous method
export repeat_evolve_until_dead, evolve_until_dead, subpop_alive

@doc """ function repeat_evolve_until_dead()
  Run repeated trails of evolve_until_dead()
  parparamd is the per-trial dictionary of paramters, and resultd is the per-trial dictionary of results.
"""
function repeat_evolve_until_dead( paramd::param_type, resultd::result_type )
  println("repeat_evolve_until_dead: N: ",paramd[:N],"  num_attributes: ",paramd[:num_attributes]," num_subpops: ",paramd[:num_subpops],"  mutStddev: ",paramd[:mutStddev])
  if paramd[:num_trials] == 1
    return evolve_until_dead( paramd, resultd )
  end
  resultd_list = result_type[]
  sum_generational_lifetime = 0.0
  sum_move_update_lifetime = 0.0
  sum_gen_limit_count = 0
  (sum_propsel_loss,sum_propsel_gain,sum_horiz_loss,sum_horiz_gain) = (0,0,0,0)
  for t = 1:paramd[:num_trials]
    resultd = evolve_until_dead( paramd, deepcopy(resultd) )
    Base.push!( resultd_list, deepcopy(resultd) )   # TODO:  what is this line doing?
    sum_generational_lifetime += resultd[:generational_lifetime] 
    sum_move_update_lifetime += resultd[:move_update_lifetime] 
    #println("repeat: mul:", resultd[:move_update_lifetime])
    sum_gen_limit_count += resultd[:gen_limit_reached_count]
    (sum_propsel_loss,sum_propsel_gain,sum_horiz_loss,sum_horiz_gain) =
        (sum_propsel_loss,sum_propsel_gain,sum_horiz_loss,sum_horiz_gain) .+ 
        (resultd[:propsel_loss],resultd[:propsel_gain],resultd[:horiz_loss],resultd[:horiz_gain])
  end
  resultd[:generational_lifetime] = sum_generational_lifetime/paramd[:num_trials]
  resultd[:move_update_lifetime] = sum_move_update_lifetime/paramd[:num_trials]
  (resultd[:propsel_loss],resultd[:propsel_gain],resultd[:horiz_loss],resultd[:horiz_gain]) = 
      (sum_propsel_loss,sum_propsel_gain,sum_horiz_loss,sum_horiz_gain) ./ paramd[:num_trials]
  resultd[:gen_limit_reached_count] = sum_gen_limit_count
  return resultd
end

@doc """ function evolve_until_dead( )
paramd is the per-trial dictionary of paramters, and resultd is the per-trial dictionary of results.
A subpop is "dead" if all members have fitness minFit.
Evolve generations until all subpops are dead, either on a single generation update, or on a move optimum update
Doesn't keep track of the other statistics.
See types.jl for the definition of param_type and result_type, and for the definition of the fields of this type.
"""
function evolve_until_dead( paramd::param_type, resultd::result_type )
  #gen_lifetime = 0
  #mu_lifetime = 0
  metapop_dead_flag = false    # used to test if all metapops since the last move updata have been dead
  gens_since_last_move_update = 0  # count generations since last move update
  #println("function evolve_until_dead")  ###
  #println("num_subpops: ",paramd[:num_subpops],"  num_emigrants: ",paramd[:num_emigrants],"  horiz_sel: ",paramd[:horiz_select],"  mutStddev: ",paramd[:mutStddev],"  topology: ",paramd[:topology])
  #println("burn_in: ",paramd[:burn_in])
  int_burn_in = Int(round(paramd[:burn_in]*paramd[:N]))
  id = [0]
  ideal = fill( paramd[:ideal_init], paramd[:num_attributes] )
  subpop_size = Int(floor(paramd[:N]/paramd[:num_subpops]))
  @assert subpop_size == paramd[:subpop_size]
  #println("N: ",resultd[:N],"  num_subpops: ",resultd[:num_subpops],"subpop_size: ",subpop_size)
  #if subpop_size*paramd[:num_subpops] != paramd[:N]
  if subpop_size*paramd[:num_subpops] != paramd[:N]
    error("N must equal subpop_size*resultd[:num_subpops]")
  end
  vt = Dict{Int64,temporal_variant_type}()
  meta_pop = init_meta_pop( paramd, vt, ideal, id )
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  resultd[:generational_lifetime]=0
  resultd[:move_update_lifetime]=0
  resultd[:gen_limit_reached_count]=0
  (propsel_loss,propsel_gain,horiz_loss,horiz_gain) = (0,0,0,0)
  g = 1
  while (resultd[:generational_lifetime]==0 || resultd[:move_update_lifetime]==0 ) && g < paramd[:ngens]+int_burn_in
    gens_since_last_move_update += 1
    subpops_alive_before = subpops_alive( meta_pop, paramd, vt )
    #println("subpops_alive_before: ",subpops_alive_before)
    if paramd[:move_time_interval] > 0 && g > int_burn_in+1 && g % paramd[:move_time_interval] == 1
      #println("move update:  metapop_dead_flag: ",metapop_dead_flag,"  gens_since_last_move_update: ",gens_since_last_move_update)
      # Check whether metapop has been dead for all generations since the last move update
      if metapop_dead_flag && resultd[:move_update_lifetime] == 0 && gens_since_last_move_update == paramd[:move_time_interval] 
        #mu_lifetime = g
        resultd[:move_update_lifetime] = g
        #println(":move_update_lifetime set to ",g)
        break
      end
      move_optima( ideal, paramd[:move_range] )
      #println(" g: ",g,"  #optimum just  moved resultd[:move_update_lifetime]: ",resultd[:move_update_lifetime])
      metapop_dead_flag = true
      #println(" g: ",g,"  #optimum just  moved   means: ", fit_means( meta_pop, vt ))
      #subpop_alive_opt_move_update( sbp, meta_pop,  vt, paramd[:minFit] )
      gens_since_last_move_update = 0
    end
    mutate_meta_pop!( meta_pop, vt, ideal, id, paramd )  # will also re-evaluate fitness
    mmeans = fit_means( meta_pop, vt )
    #println("  g: ",g," after mutate means: ",mmeans)
    for  j = 1:paramd[:num_subpops]
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )  # comment out for fitness test
    end
    subpops_alive_after_propsel= subpops_alive( meta_pop, paramd, vt )
    #println("subpops_alive_after_propsel: ",subpops_alive_after_propsel)
    #println("  g: ",g," before horiz means: ",mmeans)
    horiz_transfer( meta_pop, paramd, vt, ideal, mmeans, id, g )
    subpops_alive_after_horiz= subpops_alive( meta_pop, paramd, vt )
    #println("subpops_alive_after_horiz: ",subpops_alive_after_propsel)
    (propsel_loss, propsel_gain, horiz_loss, horiz_gain) = 
        (propsel_loss, propsel_gain, horiz_loss, horiz_gain) .+ 
        opt_gained_lost_counts( subpops_alive_before, subpops_alive_after_propsel, subpops_alive_after_horiz )
    #println("opt_counts: ",(propsel_loss, propsel_gain, horiz_loss, horiz_gain))
    mmeans, vvars, stddevs, cfvars = fit_means_vars_stddevs_cfvars( meta_pop, vt )
    #println("  g: ",g," after horiz: ",mmeans)
    #println("vt: ",vt)
    #print("g:",g,": ")
    #print_meta_pop( paramd, meta_pop, vt )
    if !metapop_alive(meta_pop,vt,paramd[:minFit]) && resultd[:generational_lifetime] == 0
      #gen_lifetime = g
      resultd[:generational_lifetime] = g
      #println(":generational_lifetime set to ",g)
    end
    if metapop_alive(meta_pop,vt,paramd[:minFit]) 
      metapop_dead_flag = false
    end
    #println("alive: ",metapop_alive(meta_pop,vt,paramd[:minFit]),"  dead_flag: ",metapop_dead_flag)
    g += 1
  end  # generational loop
  (resultd[:propsel_loss],resultd[:propsel_gain],resultd[:horiz_loss],resultd[:horiz_gain]) = 
      (propsel_loss,propsel_gain,horiz_loss,horiz_gain)
  #println("resultd gain loss: ",(resultd[:propsel_loss],resultd[:propsel_gain],resultd[:horiz_loss],resultd[:horiz_gain]))
  if resultd[:move_update_lifetime]!=0
    #println("loop terminated with move update lifetime: ",resultd[:move_update_lifetime])
  end
  if resultd[:generational_lifetime] == 0     
    #gen_lifetime = g
    resultd[:generational_lifetime] = g
    resultd[:gen_limit_reached_count] += 1
  end
  if resultd[:move_update_lifetime] == 0
    #mu_lifetime = g
    resultd[:move_update_lifetime] = g
  end
  #println("generational lifetime: ",resultd[:generational_lifetime])
  #println("move update lifetime:  ",resultd[:move_update_lifetime])
  return resultd
end

@doc """ function subpop_properties_init()
   Initializes 3 lists of Bools that describe the 'alive' properties of the subpops
"""
function subpop_properties_init( num_subpops::Int64 )
  # generational_subpop_alive[j] == true if it has been true in all previous generations and some individual has fitnes greater than minFit
  generational_subpop_alive = fill(true,num_subpops)   
  # prev_subpop_alive[j] == true if current_subpop_alive[j] was true in the previous move update
  prev_subpop_alive = fill(true,num_subpops)
  current_subpop_alive = fill(false,num_subpops)
  return subpop_properties(generational_subpop_alive,prev_subpop_alive,current_subpop_alive)
end

function subpop_alive_gen_update( sbp::subpop_properties, meta_pop::PopList,  variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  for j = 1:num_subpops
    count_indivs_at_minFit = count_individuals_at_minFit( meta_pop[j], variant_table, minFit )
    if sbp.generational_subpop_alive[j] && count_indivs_at_minFit == subpop_size
      sbp.generational_subpop_alive[j] = false  # only set to false if previously true and all individuals have fitness <= minFit
    end
    if !sbp.current_subpop_alive[j] && count_indivs_at_minFit < subpop_size
      sbp.current_subpop_alive[j] = true  # set to true if previously false and some individual has fitness greater than minFit
    end
  end
end

function subpop_alive_opt_move_update( sbp::subpop_properties, meta_pop::PopList,  variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
  if !any(sbp.prev_subpop_alive) 
    return
  end
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  for j = 1:num_subpops
    if sbp.prev_subpop_alive[j] && !sbp.current_subpop_alive[j]
      sbp.prev_subpop_alive[j] = false  # only set to false if previously true and all individuals have fitness equal to minFit
    end
    sbp.current_subpop_alive[j] = false
  end
end

@doc """ function update_subpop_alive!()    # never called
  A subpop is alive if it satisfies two conditions: it was alive in the previous generation, and it has at least one individual whose
     fitness is strictly greater than minFit.
"""
function update_subpop_alive!( subpop_alive::Vector{Bool}, prev_subpop_alive::Vector{Bool}, meta_pop::PopList, variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  for j = 1:num_subpops
    if subpop_alive[j] && count_individuals_at_minFit( meta_pop[j], variant_table, minFit ) == subpop_size
      subpop_alive[j] = false  # only set to false if previously true and all individuals have fitness less than minFit
    end
  end
end

# A subpop is alive if some individual has fitness above minFit
function subpop_alive(subpop::Population, variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
  return count_individuals_at_minFit(subpop, variant_table, minFit) != length(subpop)
end

# The metapop is alive if at least one of its subpops is alive
function metapop_alive(meta_pop::PopList, variant_table::Dict{Int64,temporal_variant_type}, minFit::Float64 )
  return count_subpops_at_minFit(meta_pop, variant_table, minFit) != length(meta_pop)
end 
