# Spatial structure simulation with horizontal transfer
export spatial_simulation, fitness

empty_variant = variant_type(-1,0.0,0,Vector{Float64}())
#vtbl = Dict{Int64,variant_type}()

@doc """ function spatial_simulation()
  Wright-Fisher model simulation (as opposed to Moran model)
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m
    mu    innovation probability
    ngens number of generations
    num_emmigrants   number of emmigrants in horizontal transfer
    num_attributes   number of quantitative attributes of a variant
    copy_dfe   Distribution of selection coefficients for the reduction in fitness during copy
    innov_dfe  Distribution of selection coefficients for innovations
    horiz_dfe  Distribution of selection coefficients for the reduction in fitness during horiz transfer
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
    quantitative==true means individuals have quantitative attributes, fitness computed by distance from ideal
    forward==true  means that horizontal transfer is done in a forward circular fashion
    extreme==true  means that horizontal transfer is done in a forward circular fashion
    neg_select==true  means that reverse proportional selection is used to select individuals to delete in horiz trans
"""
function spatial_simulation( sr::SpatialEvolution.spatial_result_type )
  variant_table = Dict{Int64,variant_type}()
  #println("spatial simulation: quantitative: ",quantitative)
  #println("sim circular_variation: ",sr.circular_variation,"  extreme_variation: ",sr.extreme_variation)
  fitness_locations = initialize_fitness_locations(sr)
  #println("fitness_locations: ",fitness_locations)
  #int_burn_in = Int(round(sr.burn_in*N+50.0))
  int_burn_in = Int(round(sr.burn_in*sr.N))   # reduce for testing
  id = Int[1]
  n = Int(floor(sr.N/sr.num_subpops))    # size of subpopulations
  println("N: ",sr.N,"  num_subpops: ",sr.num_subpops,"  n: ",n,"  use fit locations: ",sr.use_fit_locations,"  num_attributes: ",sr.num_attributes,"  ngens: ",sr.ngens,"  ne: ",sr.ne)
  cumm_means = zeros(Float64,sr.num_subpops)
  cumm_variances = zeros(Float64,sr.num_subpops)
  cumm_attr_vars = zeros(Float64,sr.num_subpops)
  pop_list = Vector{PopList}()
  subpops = PopList()
  for j = 1:sr.num_subpops
    Base.push!( subpops, Population() )
    for i = 1:n
      fit_loc_ind = fit_loc_index(sr.N,sr.num_subpops,sr.num_fit_locations,j,i)
      #println("j: ",j,"  i: ",i,"  fit_loc_ind: ",fit_loc_ind)
      Base.push!( subpops[j], new_innovation( id, 
          fit_loc_ind, sr.num_attributes, variant_table, fitness_locations ) )
    end
    #println("subpops[",j,"]: ",subpops[j] )
  end
  previous_variant_id = 1
  current_variant_id = id[1]
  Base.push!(pop_list,deepcopy(subpops))
  for g = 2:sr.ngens+int_burn_in
    previous_previous_variant_id = previous_variant_id
    previous_variant_id = current_variant_id
    current_variant_id = id[1]
    #println("g: ",g)
    for j = 1:sr.num_subpops
      for i = 1:n
        fit_loc_ind = fit_loc_index(sr.N,sr.num_subpops,sr.num_fit_locations,j,i)
        cp = copy_parent( pop_list[g-1][j][i], id, fit_loc_ind, sr.mu, 
          sr.normal_stddev, variant_table, fitness_locations )
        #println("j: ",j,"  i: ",i,"  pl: ",pop_list[g-1][j][i],"  cp: ",cp)
        subpops[j][i] = cp
      end
      #println("g:",g," j:",j,"  ",[(v,variant_table[v].attributes) for v in subpops[j]])
      subpops[j] = propsel( subpops[j], n, variant_table )
    end
    if g%2==0
      horiz_transfer_circular!( sr.N, sr.num_subpops, sr.ne, subpops, id, variant_table, fitness_locations,
          forward=true, neg_select=sr.horiz_select, emmigrant_select=sr.horiz_select )
    else
      horiz_transfer_circular!( sr.N, sr.num_subpops, sr.ne, subpops, id, variant_table, fitness_locations,
          forward=false, neg_select=sr.horiz_select, emmigrant_select=sr.horiz_select )
    end
    Base.push!(pop_list,deepcopy(subpops))
    #print_pop(STDOUT,subpops,variant_table)
    if g > int_burn_in
      (mmeans,vvars) = means(subpops,variant_table)
      cumm_means += mmeans
      cumm_variances += vvars
      avars = attr_vars(subpops,variant_table)
      cumm_attr_vars += avars
      #println("cumm_means: ",cumm_means)
      #println("cumm_variances: ",cumm_variances)
    end
    clean_up_variant_table(previous_previous_variant_id,previous_variant_id,variant_table)
  end
  cumm_means /= sr.ngens
  cumm_variances /= sr.ngens
  cumm_attr_vars /= sr.ngens
  #println("fitness mean: ",mean(cumm_means),"  variance: ",mean(cumm_variances),"  attribute_variance: ",mean(cumm_attr_vars))
  sr.fitness_mean = mean(cumm_means)
  sr.fitness_variance = mean(cumm_variances)
  sr.attribute_variance = mean(cumm_attr_vars)
  return sr
end

function fitness( attributes::Vector{Float64}, ideal::Vector{Float64} )
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  sum = 0.0
  for k = 1:length(attributes)
    sum += abs( attributes[k] - ideal[k] )
  end
  #println("fitness: attributes: ",attributes,"  ideal: ",ideal," fit: ",1.0-sum/length(attributes))
  return 1.0-sum/length(attributes)
end

function new_innovation( id::Vector{Int64}, 
    fit_loc_ind::Int64, num_attributes::Int64,
    variant_table::Dict{Int64,variant_type}, fitness_locations; quantitative::Bool=true )
  i = id[1]
  if quantitative
    #println("new innovation i: ",i,"  fit_loc_ind:",fit_loc_ind,"  num_attributes: ",num_attributes )
    variant_table[i] = variant_type( i, 0.0, fit_loc_ind, rand(num_attributes) )
    #println("variant_table[i]: ",variant_table[i])
    variant_table[i].fitness = fitness( variant_table[i].attributes, fitness_locations[fit_loc_ind].ideal )  
  end
  id[1] += 1
  i
end

@doc """  copy_parent()
  Note that the function cdfe produces an incremental change in fitness rather than a new fitness.
"""
function copy_parent( v::Int64, id::Vector{Int64}, fit_loc_ind::Int64, mu::Float64,
    normal_stddev::Float64, variant_table::Dict{Int64,SpatialEvolution.variant_type},
    fitness_locations::Vector{SpatialEvolution.fitness_location_type}; quantitative::Bool=true )
  i = id[1]
  vt = variant_table[v]
  if quantitative  # spatial structure by deviation of attributes from ideal
    vt.attributes = mutate_attributes( vt.attributes, normal_stddev )
    if rand() < mu
      innovate_attribute( vt.attributes, fit_loc_ind, fitness_locations )
    end
    #println("copy_parent v: ",v,"  fit_loc_ind: ",fit_loc_ind)
    #println("copy_parent v: ",v,"  attributes: ",vt.attributes)
    new_fit = fitness( vt.attributes, fitness_locations[fit_loc_ind].ideal )
    #=
    distance = 0.0
    println("length(vt.attributes): ",length(vt.attributes),"  length(ideal): ",length(fitness_locations[fit_loc_ind].ideal))
    for i = 1:length(vt.attributes)
      distance += abs(vt.attributes[i] - fitness_locations[fit_loc_ind].ideal[i])
    end
    new_fit = 1.0-distance/length(vt.attributes)
    =#
  else   # spatial structure by fitness increment
    ffit = cdfe()
    new_fit = vt.fitness + ffit
  end
  #println("copy_parent i: ",i,"  quantitative: ",quantitative,"  new_fit: ",new_fit)
  variant_table[i] = deepcopy(vt)
  variant_table[i].fitness = new_fit
  #variant_table[i] = variant_type(v,new_fit,vt.fitness_location,vt.attributes)  # needs to be fixed
  #println("v: ",v,"  i: ",i,"  new_fit: ",new_fit,"  vtbl[i]: ",variant_table[i].fitness)
  id[1] += 1
  return i
end  

function mutate_attributes( attributes::Vector{Float64}, normal_stddev::Float64 )
  #stddev = normal_stddev()   # Standard deviation of mutational perturbations
  #println("mutate attributes  normal_stddev: ",normal_stddev)
  #attributes = min(1.0,max(0.0,attributes+normal_stddev*randn(length(attributes))))
  for i = 1:length(attributes)
    #println("B attributes[",i,"]: ",attributes[i])
    attributes[i] += +normal_stddev*randn()
    if attributes[i] < 0
        attributes[i] += 1.0
        #println("wrapped up: ",attributes[i])
    end
    if attributes[i] > 1.0
        attributes[i] -= 1.0
        #println("wrapped down: ",attributes[i])
    end
    attributes[i] = min(1.0,max(0.0,attributes[i]))
    #println("A attributes[",i,"]: ",attributes[i])
  end
  #println("attributes: ",attributes)
  return attributes
end

function innovate_attribute( attributes::Vector{Float64}, subpop_index::Int64, fitness_locations::Vector{SpatialEvolution.fitness_location_type} )
  j = rand(1:length(attributes))   # Choose a random attribute
  #println("j: ",j,"  attribute: ",attributes[j],"  ideal: ",fitness_locations[subpop_index].ideal[j])
  attributes[j] += rand()*abs(attributes[j] - fitness_locations[subpop_index].ideal[j])*(fitness_locations[subpop_index].ideal[j]-attributes[j])
  #println("j: ",j,"  attribute: ",attributes[j],"  ideal: ",fitness_locations[subpop_index].ideal[j])
end 

function initialize_fitness_locations( sr::SpatialEvolution.spatial_result_type )
  fitness_locations = [ fitness_location_type( zeros(Float64,sr.num_attributes) ) for j = 1:sr.num_fit_locations ]
  if !sr.circular_variation && !sr.extreme_variation  # random variation---no relationship to subpop number
    #println("init sr.circular_variation: ",sr.circular_variation,"  sr.extreme_variation: ",sr.extreme_variation)
    for j = 1:sr.num_fit_locations
      for k = 1:sr.num_attributes
        if sr.ideal_min != sr.ideal_max
          fitness_locations[j].ideal[k] = sr.ideal_min+rand()*(sr.ideal_max-sr.ideal_min)
        else
          fitness_locations[j].ideal[k] = sr.ideal_min
        end
      end
    end
  elseif sr.circular_variation && !sr.extreme_variation
    increment = 2.0*(sr.ideal_max-sr.ideal_min)/sr.num_fit_locations
    mid = Int(floor(sr.num_fit_locations/2))
    for j = 1:mid
      for k = 1:sr.num_attributes
        fitness_locations[j].ideal[k] = sr.ideal_min+increment*(j-1)+(rand()*sr.ideal_range-0.5*sr.ideal_range)
      end
    end
    for j = (mid+1):sr.num_fit_locations
      for k = 1:sr.num_attributes
        fitness_locations[j].ideal[k] = sr.ideal_max-increment*(j-mid-1)+(rand()*sr.ideal_range-0.5*sr.ideal_range)
      end
    end
  elseif !sr.circular_variation && sr.extreme_variation  # randomly choose between low_value and high_value
    for j = 1:sr.num_fit_locations
      for k = 1:sr.num_attributes
        if rand() < 0.5
          fitness_locations[j].ideal[k] = sr.ideal_min+(rand()*sr.ideal_range-0.5*sr.ideal_range)
        else 
          fitness_locations[j].ideal[k] = sr.ideal_max+(rand()*sr.ideal_range-0.5*sr.ideal_range)
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
"""
function horiz_transfer_circular!( N::Int64, m::Int64, num_emmigrants::Int64, subpops::PopList, id::Vector{Int64}, 
    variant_table::Dict{Int64,variant_type}, fitness_locations::Vector{SpatialEvolution.fitness_location_type};
     forward::Bool=true, neg_select::Bool=true, emmigrant_select::Bool=true )
  n = Int(floor(N/m))    # size of subpopulations
  num_attributes = length(variant_table[subpops[1][1]].attributes)
  num_fit_locations = length(fitness_locations)
  #println("ht num_fit_locations: ",num_fit_locations)
  #println("horiz_transfer_circular! forward: ",forward,"  num_attributes: ",num_attributes)
  emmigrants = PopList()
  for j = 1:m
    if emmigrant_select
      Base.push!( emmigrants, propsel( subpops[j], num_emmigrants, variant_table ) )
    else
      Base.push!( emmigrants, subpops[j][1:num_emmigrants] )   # Neutral
    end
  end
  new_emmigrants = Population[ Population() for j = 1:m ]
  for j = 1:m
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1)
    if forward
      k = (j+m-2)%m+1
    else
      k = j%m+1
    end
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1,"  k: ",k)
    # Create new variants for the emmigrants in the new subpop
    for e in emmigrants[k]   # subpop k is the source, subpop j is the destination
      i = id[1]
      #println("e: ",e,"  i: ",i)
      #println("new emmigrant i: ",i,"  subpop_index:",k,"  num_attributes: ",num_attributes )
      variant_table[i] = deepcopy(variant_table[e])
      ii = rand(1:n)  # Choose a random index within the subpopulation
      fit_loc_ind = fit_loc_index(N,m,num_fit_locations,j,ii)
      variant_table[i].fitness_location = fit_loc_ind   # set the new fitness location
      variant_table[i].fitness = fitness( variant_table[i].attributes, fitness_locations[fit_loc_ind].ideal )  
      #variant_table[i] = variant_type( i, 0.0, j, emmigrants[j].attributes  )
      #println("variant_table[",e,"]: ",variant_table[e])
      #println("variant_table[",i,"]: ",variant_table[i])
      Base.push!( new_emmigrants[j], i )
      id[1] += 1
    end
  end
  for j = 1:m
    pop_after_deletion = Population[]
    #println("j: ",j,"  j%m+1: ",j%m+1,"  (j+m-2)%m+1: ",(j+m-2)%m+1)
    if neg_select  # use reverse proportional selection to delete elements by negative fitness
      pop_after_deletion = reverse_propsel(subpops[j],num_emmigrants,variant_table)
    else  # delete random elements to delete
      pop_after_deletion = subpops[j][1:(n-num_emmigrants)]
    end
    subpops[j] = append!( pop_after_deletion, new_emmigrants[j] )
  end
  emmigrants  # perhaps should be the modified subpops
end

function print_subpop( subpop::Vector{Int64}, variant_table::Dict{Int64,variant_type} )
  "[ $(join([ @sprintf(" %5d:%5.4f",vt,variant_table[vt].fitness) for vt in subpop ]))]"
end

function print_pop( stream::IO, subpops::PopList, variant_table::Dict{Int64,variant_type} )
  for sp in subpops
    print(stream,print_subpop(sp,variant_table))
  end
  println(stream)
end

using DataFrames
# compute and save statistics about subpopulations and populations

function means( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function attr_vars( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  num_attributes = length(variant_table[1].attributes)
  #println("attr_vars: num_attributes: ",num_attributes)
  #=
  for s in subpops
    println("subpop: ",s)
    println(s," fitness: ",[variant_table[v].fitness for v in s ],"  variance: ",var([variant_table[v].fitness for v in s ]))
    for j = 1:num_attributes
      println("attribute[",j,"]: ",[(v,variant_table[v].attributes[j]) for v in s],"  variance: ",var([variant_table[v].attributes[j] for v in s]) )
    end
  end
  =#
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

function fit_loc_index(N,num_subpops,num_fit_locs,j,i)
  n = Int(ceil(N/num_subpops))
  mult = Int(ceil(num_fit_locs/num_subpops))
  div = Int(ceil(n*num_subpops/num_fit_locs))
  return mult*(j-1) + Int(floor((i-1)/div))+1
end

function clean_up_variant_table( previous_variant_id::Int64, previous_previous_variant_id::Int64,
    variant_table::Dict{Int64,variant_type} )
  for v = previous_previous_variant_id:previous_variant_id-1
    delete!(variant_table,v)
  end
end
 
function init()
  global vtbl = Dict{Int64,variant_type}()
  global fitness_locations = fitness_location_type[ fitness_location_type(0.1,zeros(Float64,m)) ]
  #global idfe() = 1.0+rand(Distributions.Gamma(1.0,0.1))
  #global cdfe() = -rand(Distributions.Gamma(0.2,0.001))
  global N = 20
  global m = 4
  global mu = 0.2
  global cperr = 0.2
  global ne = 3  # num_emmigrants 
  global ngens = 4
end  
# Sample call to main function
# spatial_simulation( N, m, mu, cperr, ngens, ne, cdfe, idfe, vtbl, fitness_locations )
