export repeat_evolve_until_dead, evolve_until_dead

@doc """ function repeat_evolve_until_dead()
  Run repeated trails of evolve_until_dead()
"""
function repeat_evolve_until_dead( tr::temporal_result_type )
  println("repeat_evolve_until_dead: num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutStddev: ",tr.mutStddev,
        "  topology: ",tr.topology,"  linfit_slope: ",tr.linfit_slope)   ###
  if tr.num_trials == 1
    return evolve( tr )
  end
  tr_list = temporal_result_type[]
  sum_generational_lifetime = 0.0
  sum_move_update_lifetime = 0.0
  sum_gen_limit_count = 0
  for t = 1:tr.num_trials
    tr = evolve_until_dead( tr )
    Base.push!( tr_list, deepcopy(tr) )
    sum_generational_lifetime += tr.generational_lifetime 
    sum_move_update_lifetime += tr.move_update_lifetime 
    sum_gen_limit_count += tr.gen_limit_reached_count
  end
  tr.generational_lifetime = sum_generational_lifetime/tr.num_trials
  tr.move_update_lifetime = sum_move_update_lifetime/tr.num_trials
  tr.gen_limit_reached_count = sum_gen_limit_count
  return tr
ed

end
@doc """ function evolve_until_dead( )
A subpop is "dead" if all members have fitness minFit.
Evolve generations until all subpops are dead, either on a single generation update, or on a move optimum update
Doesn't keep track of the other statistics.
See types.jl for the definition of temporal_result_type, and for the definition of the fields of this type.
"""
function evolve_until_dead( tr::temporal_result_type )
  println("function evolve_until_dead")  ###
  #println("num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutStddev: ",tr.mutStddev,"  topology: ",tr.topology)
  int_burn_in = Int(round(tr.burn_in*tr.N))
  id = [0]
  ideal = fill( tr.ideal_init, tr.num_attributes )
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  if subpop_size*tr.num_subpops != tr.N
    error("N must equal subpop_size*tr.num_subpops")
  end
  vt = Dict{Int64,variant_type}()
  meta_pop = init_meta_pop( tr, vt, ideal, id )
  #subpop_alive = fill(true,tr.num_subpops)  # Should be used with uniform start
  #prev_subpop_alive = fill(true,tr.num_subpops)  # Should be used with uniform start
  sbp = subpop_properties_init( tr.num_subpops)
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  tr.generational_lifetime=0
  tr.move_update_lifetime=0
  tr.gen_limit_reached_count=0
  g = 1
  #while (any(sbp.generational_subpop_alive) || any(sbp.prev_subpop_alive)) && g < tr.ngens+int_burn_in
  while (tr.generational_lifetime==0 || tr.move_update_lifetime==0 ) && g < tr.ngens+int_burn_in
    if tr.move_time_interval > 0 && g > int_burn_in && g % tr.move_time_interval == 0
      #println(" g: ",g,"  #optimum to be moved   means: ", fmeans( meta_pop, vt ))
      #println("B g: ",g,"  sbp.prev_subpop_alive: ",sbp.prev_subpop_alive )
      #println("B g: ",g,"  sbp.current_subpop_alive: ",sbp.current_subpop_alive )
      move_optima( ideal, tr.move_range )
      #println(" g: ",g,"  #optimum just  moved   means: ", fmeans( meta_pop, vt ))
      subpop_alive_opt_move_update( sbp, meta_pop,  vt, tr.minFit )
      #println("A g: ",g,"  sbp.prev_subpop_alive: ",sbp.prev_subpop_alive )
      #println("A g: ",g,"  sbp.current_subpop_alive: ",sbp.current_subpop_alive )
    end
    mutate_meta_pop!( meta_pop, vt, ideal, id, tr )  # will also re-evaluate fitness
    for  j = 1:tr.num_subpops
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )  # comment out for fitness test
    end
    mmeans = fmeans( meta_pop, vt )
    #println("  g: ",g," before horiz: ",mmeans)
    horiz_transfer( meta_pop, tr, vt, ideal, mmeans, id, g )
    mmeans, vvars = means_vars( meta_pop, vt )
    #println("  g: ",g," after horiz: ",mmeans)
    #println("vt: ",vt)
    subpop_alive_gen_update( sbp, meta_pop, vt, tr.minFit )
    #println("  g: ",g,"  sbp.generational_subpop_alive: ",sbp.generational_subpop_alive )
    #println("  g: ",g,"  sbp.prev_subpop_alive: ",sbp.prev_subpop_alive )
    #println("  g: ",g,"  sbp.current_subpop_alive: ",sbp.current_subpop_alive )
    if tr.generational_lifetime==0 && !any(sbp.generational_subpop_alive)
      #println("g: ",g," All subpops died (single generation update).")
      tr.generational_lifetime = g
    elseif tr.move_update_lifetime==0 && !any(sbp.prev_subpop_alive)
      #println("g: ",g," All subpops died (move optimum update).")
      tr.move_update_lifetime = g
    end
    g += 1
  end
  if g == tr.ngens+int_burn_in
    #println("evolve_until_dead ran for ",g," generations when the generation limit was reached")
    tr.gen_limit_reached_count += 1
    if tr.generational_lifetime > 0
      #println("all subpops died at generation ",tr.generational_lifetime,"  (generational update).")
    else
      tr.generational_lifetime = g
    end
    if tr.move_update_lifetime > 0
      #println("all subpops died at generation ",tr.move_update_lifetime,"(move optimum update).")
    else
      tr.move_update_lifetime = g
    end
  else 
    #println("all subpops died at generation ",tr.generational_lifetime,"  (genupdate).")
    #println("all subpops died at generation ",tr.move_update_lifetime,"  (move optimum update).")
  end
  return tr
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

