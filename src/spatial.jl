# Spatial structure simulation with horizontal transfer
export spatial_simulation, repeat_spatial, fitness, print_fit_locations

@doc """ function repeat_spatial( )
Runs spatial_simulation()  sr[:num_subpops] times, and averages the results.
See types.jl for the definition of param_type and result_type, and see types.jl for legal keys for this type.
"""
function repeat_spatial( sp::param_type, sr::result_type )  # sp is parameter dictionary, sr is results dictionary
  #=
  if sp[:num_trials] == 1
    return spatial_simulation( sp, sr )
  end
  =#
  #int_burn_in = Int(round(sp[:burn_in]*sp[:N]))
  #println("int_burn_in: ",int_burn_in)
  sr_list = result_type[]
  sum_fit_mean = 0.0
  sum_fit_vars = 0.0
  sum_fit_stddev = 0.0
  sum_fit_coef_var = 0.0
  sum_attr_mean = 0.0
  sum_attr_vars = 0.0
  sum_attr_stddev = 0.0
  sum_attr_coef_var = 0.0

  for t = 1:sp[:num_trials]
    sr = spatial_simulation( sp, sr ) 
    Base.push!( sr_list, deepcopy(sr) )
    sum_fit_mean += sr[:fitness_mean]
    sum_fit_vars += sr[:fitness_variance]
    sum_fit_stddev += sr[:fitness_stddev]
    sum_fit_coef_var += sr[:fitness_coef_var]
    sum_attr_mean += sr[:attribute_mean]
    sum_attr_vars += sr[:attribute_variance]
    sum_attr_stddev += sr[:attribute_stddev]
    sum_attr_coef_var += sr[:attribute_coef_var]
  end

  sr[:fitness_mean] = sum_fit_mean/sp[:num_trials]
  sr[:fitness_variance] = sum_fit_vars/sp[:num_trials]
  sr[:fitness_stddev] = sum_fit_stddev/sp[:num_trials]
  sr[:fitness_coef_var] = sum_fit_coef_var/sp[:num_trials]
  sr[:attribute_mean] = sum_attr_vars/sp[:num_trials]
  sr[:attribute_variance] = sum_attr_vars/sp[:num_trials]
  sr[:attribute_stddev] = sum_attr_stddev/sp[:num_trials]
  sr[:attribute_coef_var] = sum_attr_coef_var/sp[:num_trials]
  return sr
end

empty_variant = temporal_variant_type(0.0,Vector{Float64}(),0)

@doc """ function spatial_simulation()
  Wright-Fisher model simulation (as opposed to Moran model)
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m
    mu    innovation probability
    ngens number of generations
    num_emigrants   number of emigrants in horizontal transfer
    num_attributes   number of quantitative attributes of a variant
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
    quantitative==true means individuals have quantitative attributes, fitness computed by distance from ideal
    forward==true  means that horizontal transfer is done in a forward circular fashion
    extreme==true  means that horizontal transfer is done in a forward circular fashion
    neg_select==true  means that reverse proportional selection is used to select individuals to delete in horiz trans
"""
function spatial_simulation( sp::TemporalEvolution.param_type, sr::TemporalEvolution.result_type )
  vt = Dict{Int64,temporal_variant_type}()
  #println("sim circular_variation: ",sp[:circular_variation],"  extreme_variation: ",sp[:extreme_variation])
  fitness_locations = initialize_fitness_locations(sp)
  #print_fit_locations( fitness_locations )
  #println("fitness_locations: ",fitness_locations)
  if typeof(sp[:burn_in]) == Int64
    int_burn_in = sp[:burn_in] 
  else
    int_burn_in = Int(round(sp[:burn_in]*sp[:N]))  # moved to run_spatial.jl 
  end
  #println("int_burn_in: ",sp[:int_burn_in])
  id = Int[1]
  n = Int(floor(sp[:N]/sp[:num_subpops]))    # size of subpopulations
  #println("N: ",sp[:N],"  num_subpops: ",sp[:num_subpops],"  n: ",n,"  use fit locations: ",sp[:use_fit_locations],"  num_attributes: ",sp[:num_attributes],"  ngens: ",sp[:ngens],"  ne: ",sp[:num_emigrants])
  cumm_fit_means = zeros(sp[:num_subpops])
  cumm_fit_vars = zeros(sp[:num_subpops])
  cumm_fit_stddevs = zeros(sp[:num_subpops])
  cumm_fit_coef_vars = zeros(sp[:num_subpops])
  cumm_attr_means = zeros(sp[:num_subpops])
  cumm_attr_vars = zeros(sp[:num_subpops])
  cumm_attr_stddevs = zeros(sp[:num_subpops])
  cumm_attr_coef_vars = zeros(sp[:num_subpops])
  #pop_list = Vector{PopList}()
  subpops = PopList()
  for j = 1:sp[:num_subpops]
    Base.push!( subpops, Population() )
    for i = 1:n
      fit_loc_ind = fit_loc_index(sp[:N],sp[:num_subpops],sp[:num_fit_locations],j,i)
      #println("j: ",j,"  i: ",i,"  fit_loc_ind: ",fit_loc_ind)
      Base.push!( subpops[j], new_innovation( id, 
          fit_loc_ind, sp[:num_attributes], sp[:neutral], sp[:linfit_slope], vt, fitness_locations ) )
    end
    #println("subpops[",j,"]: ",subpops[j], "  pop attr: ",[ vt[subpops[j][i]].attributes[1] for i = 1:n ] )
  end
  previous_variant_id = 1
  current_variant_id = id[1]
  #Base.push!(pop_list,deepcopy(subpops))
  previous_subpops = deepcopy(subpops)
  count_gens = 0
  for g = 1:sp[:ngens]+int_burn_in
    #println("before g: ",g,"  pop: ",subpops[1],"  pop attr: ",[ vt[subpops[1][i]].attributes[1] for i = 1:n ])
    previous_previous_variant_id = previous_variant_id
    previous_variant_id = current_variant_id
    current_variant_id = id[1]
    #println("g: ",g,"  int_burn_in: ",int_burn_in,"  count_gens: ",count_gens)
    for j = 1:sp[:num_subpops]
      for i = 1:n
        fit_loc_ind = fit_loc_index(sp[:N],sp[:num_subpops],sp[:num_fit_locations],j,i)
        cp = copy_parent( previous_subpops[j][i], id, fit_loc_ind, sp[:mu], sp[:mutStddev], sp[:additive_error], sp[:neutral], sp[:linfit_slope],
            vt, fitness_locations )
        subpops[j][i] = cp
      end
      #println("g:",g," j:",j,"  ",[(v,vt[v].attributes) for v in subpops[j]])
      if !sp[:neutral]
        subpops[j] = propsel( subpops[j], n, vt )
      else # Wright-Fisher copy
        r = rand(1:n,n)
        subpops[j] = subpops[j][ r ]
      end
      #println("g: ",g,"  pop: ",subpops[j],"  pop attr: ",[ vt[subpops[j][i]].attributes[1] for i = 1:n ],
      #    "  pop fitnesses: ",[ vt[subpops[j][i]].fitness for i = 1:n ])
    end
    if sp[:num_emigrants] > 0 && g%2==0
      horiz_transfer_circular!( sp[:N], sp[:num_subpops], sp[:num_emigrants], subpops, id, vt, fitness_locations, sp[:linfit_slope],
          forward=true, neg_select=sp[:horiz_select], emigrant_select=sp[:horiz_select], neutral=sp[:neutral] )
    elseif sp[:num_emigrants] > 0
      horiz_transfer_circular!( sp[:N], sp[:num_subpops], sp[:num_emigrants], subpops, id, vt, fitness_locations, sp[:linfit_slope],
          forward=false, neg_select=sp[:horiz_select], emigrant_select=sp[:horiz_select], neutral=sp[:neutral] )
    end
    previous_subpops = deepcopy(subpops)
    #print_pop(STDOUT,subpops,vt)
    if g > int_burn_in
      count_gens += 1
      mmeans, vvars, stddevs, cfvars = fit_means_vars_stddevs_cfvars( subpops, vt )
      cumm_fit_means += mmeans
      cumm_fit_vars += vvars
      cumm_fit_stddevs += stddevs
      cumm_fit_coef_vars += cfvars
      mmeans, vvars, stddevs, cfvars = attr_means_vars_stddevs_cfvars( subpops, vt )
      cumm_attr_means += mmeans
      cumm_attr_vars += vvars
      cumm_attr_stddevs += stddevs
      cumm_attr_coef_vars += cfvars
    end
    clean_up_variant_table(previous_previous_variant_id,previous_variant_id,vt)
  end  # for g
  sr[:fitness_mean] = mean(cumm_fit_means/sp[:ngens])
  sr[:fitness_variance] = mean(cumm_fit_vars/sp[:ngens])
  sr[:fitness_stddev] = mean(cumm_fit_stddevs/sp[:ngens])
  sr[:fitness_coef_var] = mean(cumm_fit_coef_vars/sp[:ngens])
  sr[:attribute_mean] = mean(cumm_attr_means/sp[:ngens])
  sr[:attribute_variance] = mean(cumm_attr_vars/sp[:ngens])
  sr[:attribute_stddev] = mean(cumm_attr_stddevs/sp[:ngens])
  sr[:attribute_coef_var] = mean(cumm_attr_coef_vars/sp[:ngens])
  #println("sr[:fitness_mean] ",sr[:fitness_mean])
  #println("sr[:fitness_variance] ",sr[:fitness_variance])
  #println("sr[:attribute_variance] ",sr[:attribute_variance])
  #println("sr[:attribute_coef_var] ",sr[:attribute_coef_var])

  return sr
end

function fitness( attributes::Vector{Float64}, ideal::Vector{Float64}, neutral::Bool, linfit_slope::Float64 )
  if neutral
    return 1.0
  end
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  #println("fitness: attr: ",attributes)
  dis = 0.0
  for k = 1:length(attributes)
    dis += abs( attributes[k] - ideal[k] )
  end
  #println("fitness: attributes: ",attributes,"  ideal: ",ideal," fit: ",1.0-dis/length(attributes))
  if linfit_slope == 0.0  # use the older linear method of computing fitness
    result = max(1.0-dis/length(attributes),0.0)
  else  # inverse method of computing fitness added on 11/14/17
    result = 1.0/(linfit_slope*dis+1.0)
  end
  #println("fitness: attr: ",attributes,"  result: ",result)
  return result
end

function new_innovation( id::Vector{Int64}, 
    fit_loc_ind::Int64, num_attributes::Int64, neutral::Bool, linfit_slope::Float64,
    variant_table::Dict{Int64,temporal_variant_type}, fitness_locations; quantitative::Bool=true )
  #println("new innov: fit_loc_ind: ",fit_loc_ind)
  #println("new innov: ideal: ",fitness_locations[fit_loc_ind].ideal)
  i = id[1]
  if quantitative
    #println("new innovation i: ",i,"  fit_loc_ind:",fit_loc_ind,"  num_attributes: ",num_attributes )
    variant_table[i] = temporal_variant_type( 0.0, fitness_locations[fit_loc_ind].ideal, fit_loc_ind )
    #println("variant_table[i]: ",variant_table[i])
    variant_table[i].fitness = fitness( variant_table[i].attributes, fitness_locations[fit_loc_ind].ideal, neutral, linfit_slope )  
    #println("new innov: fitness: ",variant_table[i].fitness)
  end
  id[1] += 1
  i
end

@doc """  copy_parent()
"""
function copy_parent( v::Int64, id::Vector{Int64}, fit_loc_ind::Int64, mu::Float64,
    mutStddev::Float64, additive_error::Bool, neutral::Bool, linfit_slope::Float64, variant_table::Dict{Int64,TemporalEvolution.temporal_variant_type},
    fitness_locations::Vector{TemporalEvolution.fitness_location_type} )
  i = id[1]
  vt = variant_table[v]
  new_attributes = mutate_attributes( vt.attributes, mutStddev, additive_error )
  if mu > 0 && rand() < mu
    innovate_attribute( new_attributes, fit_loc_ind, fitness_locations )
  end
  #println("copy_parent v: ",v,"  fit_loc_ind: ",fit_loc_ind)
  #println("copy_parent v: ",v,"  attributes: ",vt.attributes,"  new_attr: ",new_attributes)
  new_fit = fitness( new_attributes, fitness_locations[fit_loc_ind].ideal, neutral, linfit_slope )
  #println("copy_parent i: ",i,"  quantitative: ",quantitative,"  new_fit: ",new_fit)
  variant_table[i] = deepcopy(vt)
  variant_table[i].fitness = new_fit
  variant_table[i].attributes = new_attributes
  #variant_table[i] = temporal_variant_type(v,new_fit,vt.fitness_location,vt.attributes)  # needs to be fixed
  #println("v: ",v,"  i: ",i,"  new_fit: ",new_fit,"  vtbl[i]: ",variant_table[i].fitness)
  id[1] += 1
  return i
end  

function mutate_attributes( attributes::Vector{Float64}, mutStddev::Float64, additive_error::Bool )
  #stddev = mutStddev()   # Standard deviation of mutational perturbations
  #println("mutate attributes  mutStddev: ",mutStddev)
  #attributes = min(1.0,max(0.0,attributes+mutStddev*randn(length(attributes))))
  new_attributes = deepcopy(attributes)
  if additive_error
    for i = 1:length(new_attributes)
      #println("B attributes[",i,"]: ",new_attributes[i])
      new_attributes[i] += +mutStddev*randn()
      if new_attributes[i] < 0
          new_attributes[i] += 1.0
          #println("wrapped up: ",new_attributes[i])
      end
      if new_attributes[i] > 1.0
          new_attributes[i] -= 1.0
          #println("wrapped down: ",new_attributes[i])
      end
      new_attributes[i] = min(1.0,max(0.0,new_attributes[i]))
      #println("A attributes[",i,"]: ",new_attributes[i])
    end
  else  # multiplicative copy error
    for i = 1:length(new_attributes)
      if new_attributes[i] <= 0.0
        #println("neg attribute: ",new_attributes[i])
        new_attributes[i] = 1.0e-6
      end
      multiplier = (1.0+mutStddev*randn())
      #println("multiplier: ",multiplier)
      while multiplier <= 1.0e-6
        #println("neg multiplier")
        multiplier = (1.0+mutStddev*randn())
      end
      new_attributes[i] *= multiplier
      if new_attributes[i] < 0.0
        #println("negative attribute with i=",i,": attribute: ",new_attributes[i])
      end
    end
  end
  #println("mutate: new_attributes: ",new_attributes)
  return new_attributes
end

# Changed to the "linear" model of innovation from the "quadratic" model on 2/14/18.
function innovate_attribute( attributes::Vector{Float64}, subpop_index::Int64, fitness_locations::Vector{TemporalEvolution.fitness_location_type} )
  j = rand(1:length(attributes))   # Choose a random attribute
  Bdiff = abs(attributes[j]-fitness_locations[subpop_index].ideal[j])
  #println("B  j: ",j,"  attribute: ",attributes[j],"  ideal: ",fitness_locations[subpop_index].ideal[j]," diff: ",Bdiff)
  #attributes[j] += rand()*abs(attributes[j] - fitness_locations[subpop_index].ideal[j])*(fitness_locations[subpop_index].ideal[j]-attributes[j])
  attributes[j] += rand()*(fitness_locations[subpop_index].ideal[j] - attributes[j])   # additive innovation
  #=
  r = rand()
  println("r: ",r,"  ratio: ",fitness_locations[subpop_index].ideal[j]/attributes[j])
  attributes[j] *= rand()*(fitness_locations[subpop_index].ideal[j]/attributes[j])   # multiplicative innovation
  =#
  Adiff = abs(attributes[j]-fitness_locations[subpop_index].ideal[j])
  #println("A  j: ",j,"  attribute: ",attributes[j],"  ideal: ",fitness_locations[subpop_index].ideal[j]," diff: ",Adiff)
  #println("Ddiff: ",Bdiff-Adiff)
  @assert(Bdiff-Adiff>0.0)
end 

@doc """ function initialize_fitness_locations()
  Sets up ideal value for each fitness location.
  There are 3 kinds of spatial fitness variation:
  1)  Random:  ideal is chosen randomly between sp[:ideal_min] and sp[:ideal_max]
  2)  Circular:  as one iterates through fit locations, ideal starts close to sp[:ideal_min], then moves toward
      sp[:ideal_max], then moves back to sp[:ideal_min].
  3)  Extreme:  ideal is chosen to be within 0.5*sp[:ideal_range] of sp[:ideal_min] or sp[:ideal_max]
"""
function initialize_fitness_locations( sp::param_type )
  #println("initialize_fitness_locations: num_fit_locations: ",sp[:num_fit_locations] )
  fitness_locations = [ fitness_location_type( zeros(Float64,sp[:num_attributes]) ) for j = 1:sp[:num_fit_locations] ]
  #println("init sp[:circular_variation]: ",sp[:circular_variation],"  sp[:extreme_variation]: ",sp[:extreme_variation])
  if !sp[:circular_variation] && !sp[:extreme_variation]  # random variation---no relationship to subpop number
    for j = 1:sp[:num_fit_locations]
      for k = 1:sp[:num_attributes]
        if sp[:ideal_min] != sp[:ideal_max]
          fitness_locations[j].ideal[k] = sp[:ideal_min]+rand()*(sp[:ideal_max]-sp[:ideal_min])
        else
          fitness_locations[j].ideal[k] = sp[:ideal_min]
        end
      end
      #println("j: ",j,"  ideal: ",fitness_locations[j].ideal)
    end
  elseif sp[:circular_variation] && !sp[:extreme_variation]
    increment = 2.0*(sp[:ideal_max]-sp[:ideal_min])/sp[:num_fit_locations]
    mid = Int(floor(sp[:num_fit_locations]/2))
    for j = 1:mid
      for k = 1:sp[:num_attributes]
        fitness_locations[j].ideal[k] = sp[:ideal_min]+increment*(j-1)+(rand()*sp[:ideal_range]-0.5*sp[:ideal_range])
      end
    end
    for j = (mid+1):sp[:num_fit_locations]
      for k = 1:sp[:num_attributes]
        fitness_locations[j].ideal[k] = sp[:ideal_max]-increment*(j-mid-1)+(rand()*sp[:ideal_range]-0.5*sp[:ideal_range])
      end
      #println("j: ",j,"  ideal: ",fitness_locations[j].ideal)
    end
  elseif !sp[:circular_variation] && sp[:extreme_variation]  # randomly choose between low_value and high_value
    for j = 1:sp[:num_fit_locations]
      for k = 1:sp[:num_attributes]
        if rand() < 0.5
          fitness_locations[j].ideal[k] = sp[:ideal_min]+(rand()*sp[:ideal_range]-0.5*sp[:ideal_range])
        else 
          fitness_locations[j].ideal[k] = sp[:ideal_max]+(rand()*sp[:ideal_range]-0.5*sp[:ideal_range])
        end
      end
      #println("j: ",j,"  ideal: ",fitness_locations[j].ideal)
    end
  else
    error("combining circular variation and extreme variation in initialize_fitness_locations is not legal.")
  end
  return fitness_locations
end  

@doc """ horiz_transfer_circular!()
  Transfers variants between subpopulations in a circular fashion (either forward or backward).
  Elements to be transfered are selected by proportional selection.
  Elements to be replaced can be random or selected by reverse proportional selection depending on the flag neg_select.
  subpops is modified by this function (as a side effect)
  Note:  m  is the number of subpops
"""
function horiz_transfer_circular!( N::Int64, m::Int64, num_emigrants::Int64, subpops::PopList, id::Vector{Int64}, 
    variant_table::Dict{Int64,temporal_variant_type}, fitness_locations::Vector{TemporalEvolution.fitness_location_type}, linfit_slope::Float64;
     forward::Bool=true, neg_select::Bool=true, emigrant_select::Bool=true, neutral::Bool=false )
  n = Int(floor(N/m))    # size of subpopulations
  num_attributes = length(variant_table[subpops[1][1]].attributes)
  num_fit_locations = length(fitness_locations)
  #println("ht num_fit_locations: ",num_fit_locations)
  #println("horiz_transfer_circular! forward: ",forward,"  num_attributes: ",num_attributes)
  emigrants = PopList()
  for j = 1:m
    if emigrant_select
      Base.push!( emigrants, propsel( subpops[j], num_emigrants, variant_table ) )
    else
      # TODO:  do random choice instead of first elements
      # Added shuffle on 1/23/18.  Not efficient, but should be correct.
      Base.push!( emigrants, Base.shuffle(subpops[j])[1:num_emigrants] )   # Neutral
    end
  end
  new_emigrants = Population[ Population() for j = 1:m ]
  for j = 1:m
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1)
    if forward
      k = (j+m-2)%m+1
    else
      k = j%m+1
    end
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1,"  k: ",k)
    # Create new variants for the emigrants in the new subpop
    for e in emigrants[k]   # subpop k is the source, subpop j is the destination
      i = id[1]
      #println("e: ",e,"  i: ",i)
      #println("new emigrant i: ",i,"  subpop_index:",k,"  num_attributes: ",num_attributes )
      variant_table[i] = deepcopy(variant_table[e])
      ii = rand(1:n)  # Choose a random index within the subpopulation
      fit_loc_ind = fit_loc_index(N,m,num_fit_locations,j,ii)
      variant_table[i].fitness_location = fit_loc_ind   # set the new fitness location
      variant_table[i].fitness = fitness( variant_table[i].attributes, fitness_locations[fit_loc_ind].ideal, neutral, linfit_slope )  
      #variant_table[i] = temporal_variant_type( i, 0.0, j, emigrants[j].attributes  )
      #println("variant_table[",e,"]: ",variant_table[e])
      #println("variant_table[",i,"]: ",variant_table[i])
      Base.push!( new_emigrants[j], i )
      id[1] += 1
    end
  end
  for j = 1:m
    pop_after_deletion = Population[]
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1)
    if neg_select  # use reverse proportional selection to delete elements by negative fitness
      pop_after_deletion = reverse_propsel(subpops[j],num_emigrants,variant_table)
    else  # delete random elements to delete
      pop_after_deletion = subpops[j][1:(n-num_emigrants)]
    end
    subpops[j] = append!( pop_after_deletion, new_emigrants[j] )
  end
  emigrants  # perhaps should be the modified subpops
end

function print_subpop( subpop::Vector{Int64}, variant_table::Dict{Int64,temporal_variant_type} )
  "[ $(join([ @sprintf(" %5d:%5.4f",vt,variant_table[vt].fitness) for vt in subpop ]))]"
end

function print_pop( stream::IO, subpops::PopList, variant_table::Dict{Int64,temporal_variant_type} )
  for sp in subpops
    print(stream,print_subpop(sp,variant_table))
  end
  println(stream)
end

function print_fit_locations( fitness_locations::Array{TemporalEvolution.fitness_location_type,1} )
  num_fit_locations = length(fitness_locations)
  num_attributes = length(fitness_locations[1].ideal)
  println("fit_locations:  num_fit_locations: ",num_fit_locations,"  num_attributes: ",num_attributes)
  for i = 1:num_fit_locations
    #print("subpop ",i,":")
    for j = 1:num_attributes
      @printf(" %6.4f",fitness_locations[i].ideal[j])
    end
    println()
  end
end

# compute and save statistics about subpopulations and populations
function means( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function attr_vars( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type}, num_attributes::Int64 )
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

function attr_coef_vars( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type}, num_attributes::Int64 )
  #for s in subpops
  #  println("attr_vars: attributes: ",[[variant_table[v].attributes[j] for v in s] for j =1:num_attributes])
  #end
  avg_coef_vars = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_coef_vars = [ coef_var([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    #println(s," att_coef_vars: ",att_coef_vars)
    avg_coef_vars[i] = mean(att_coef_vars)
    i += 1
  end
  println("avg_coef_vars: ",avg_coef_vars)
  return avg_coef_vars
end

function fit_loc_index(N,num_subpops,num_fit_locs,j,i)
  n = Int(ceil(N/num_subpops))
  mult = Int(ceil(num_fit_locs/num_subpops))
  div = Int(ceil(n*num_subpops/num_fit_locs))
  #println("fit_loc_index num_subpops: ",num_subpops,"  num_fit_locs: ",num_fit_locs,"  j: ",j,"  i: ",i,"  result: ",mult*(j-1) + Int(floor((i-1)/div))+1)
  return mult*(j-1) + Int(floor((i-1)/div))+1
end

function clean_up_variant_table( previous_variant_id::Int64, previous_previous_variant_id::Int64,
    variant_table::Dict{Int64,temporal_variant_type} )
  #println("clean up:  ppv: ",previous_previous_variant_id,"  pv: ",previous_variant_id)
  for v = previous_variant_id:previous_previous_variant_id-1
    delete!(variant_table,v)
  end
end
