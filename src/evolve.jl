export temporal_result, run_evolve, evolve, horiz_transfer_circular!, mutate_meta_pop!
    
#include("types.jl")
#include("propsel.jl")
#include("fitness.jl") 

function temporal_result( N::Int64, num_attributes::Int64, num_subpops::Int64, ngens::Int64, mutation_stddev::Float64, num_emmigrants::Int64,
    move_range::Float64, move_time_interval::Int64, opt_loss_cutoff::Float64, horiz_select::Bool=false, uniform_start::Bool=true )
  burn_in = 1.0
  ideal_init = 0.5
  return temporal_result_type( N, num_subpops, num_emmigrants, num_attributes, ngens, burn_in, uniform_start, horiz_select, mutation_stddev,
      ideal_init, move_range, move_time_interval, opt_loss_cutoff, 0.0, 0.0, 0.0, 0.0 )
end

function run_evolve(  N::Int64, num_attributes::Int64, num_subpops::Int64, ngens::Int64, mutation_stddev::Float64, num_emmigrants::Int64,
      move_range::Float64, move_time_interval::Int64, opt_loss_cutoff::Float64, horiz_select::Bool=false, uniform_start::Bool=true )
  tr = temporal_result( N, num_attributes, num_subpops, ngens, mutation_stddev, num_emmigrants, move_range, move_time_interval, opt_loss_cutoff, 
      horiz_select, uniform_start )
  evolve( tr )
end


function evolve( tr::temporal_result_type )
  println("num_subpops: ",tr.num_subpops,"  num_emmigrants: ",tr.ne,"  horiz_sel: ",tr.horiz_select,"  mutation stddev: ",tr.mutation_stddev)
  int_burn_in = Int(round(tr.burn_in*tr.N)) + 10
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
    vt[1] = variant_type( fitness( attr, ideal ), attr )
    id[1] += 1
  else
    meta_pop = [ (j-1)*subpop_size+collect(1:subpop_size)  for j = 1:tr.num_subpops ]
    for i = 1:tr.N
      attr = rand(tr.num_attributes)
      vt[i] = variant_type( fitness( attr, ideal ), attr )
    end
    id[1] += tr.N
  end
  #println("meta_pop: ",[meta_pop[j] for j = 1:length(meta_pop)])
  #println("vt: ",vt)
  for g = 2:(tr.ngens+int_burn_in)
    if g > int_burn_in && g % tr.move_time_interval == 0
      #print("g: ",g,"  mmeans: ",mmeans)
      #print("   fstdev: ",sqrt(vvars))
      #print("   astdev: ",sqrt(att_vars))
      #println("  count below cutoff: ",count_pops_below_fit_cutoff( mmeans, tr.opt_loss_cutoff ))
      #println("  count below cutoff: ",count_below_cutoff,"   cumm_count_below_cutoff: ",cumm_count_below_cutoff)
      move_optima( ideal, tr.move_range )
    end
    mutate_meta_pop!( meta_pop, vt, ideal, id, tr.mutation_stddev )  # will also re-evaluate fitness
    for  j = 1:tr.num_subpops
      meta_pop[j] = propsel( meta_pop[j], subpop_size, vt )
    end
    if g%2==0
      horiz_transfer_circular!( meta_pop, tr, vt, ideal, id,
          forward=true, neg_select=tr.horiz_select, emmigrant_select=tr.horiz_select )
    else
      horiz_transfer_circular!( meta_pop, tr, vt, ideal, id,
          forward=false, neg_select=tr.horiz_select, emmigrant_select=tr.horiz_select )
    end
    if g > int_burn_in
      (mmeans, vvars) = means( meta_pop, vt )
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
  tr.mean_fraction_subpops_below_cutoff = cumm_count_below_cutoff/tr.num_subpops/tr.ngens
  tr.fitness_mean = mean(cumm_means/tr.ngens)
  #println("cumm_means/tr.ngens: ",cumm_means/tr.ngens,"  mean: ",tr.fitness_mean)
  #println("cumm_count_below_cutoff: ",cumm_count_below_cutoff)
  #println("cumm_count_below_cutoff/num_subpops: ",cumm_count_below_cutoff/tr.num_subpops)
  tr.fitness_variance = mean(cumm_vars/tr.ngens)
  tr.attribute_variance = mean(cumm_attr_vars/tr.ngens)
  return tr

end

@doc """ horiz_transfer_circular!()
  Transfers variants between subpopulations in a circular fashion (either forward or backward).
  Elements to be transfered are selected by proportional selection.
  Elements to be replaced can be random or selected by reverse proportional selection depending on the flag neg_select.
  subpops is modified by this function (as a side effect)
"""
function horiz_transfer_circular!(  meta_pop::PopList, tr::temporal_result_type, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64};
     forward::Bool=true, neg_select::Bool=true, emmigrant_select::Bool=true )
  #println("horiz_transfer_circular! forward: ",forward,"  num_attributes: ",tr.num_attributes)
  subpop_size = Int(floor(tr.N/tr.num_subpops))
  emmigrants = PopList()
  for j = 1:tr.num_subpops
    if emmigrant_select
      Base.push!( emmigrants, propsel( meta_pop[j], tr.ne, vt ) )
    else
      Base.push!( emmigrants, meta_pop[j][1:tr.ne] )   # Neutral
    end
  end
  #println("emmigrants: ",emmigrants)
  new_emmigrants = Population[ Population() for j = 1:tr.num_subpops ]
  for j = 1:tr.num_subpops
    #println("j: ",j,"  j%tr.num_subpops+1: ",j%tr.num_subpops+1,"  (j+tr.num_subpops-2)%tr.num_subpops+1: ",(j+tr.num_subpops-2)%tr.num_subpops+1)
    if forward
      k = (j+tr.num_subpops-2)%tr.num_subpops+1
    else
      k = j%tr.num_subpops+1
    end
    #println("j: ",j,"  j%tr.num_subpops+1: ",j%tr.num_subpops+1,"  (j+tr.num_subpops-2)%tr.num_subpops+1: ",(j+tr.num_subpops-2)%tr.num_subpops+1,"  k: ",k)
    # Create new variants for the emmigrants in the new subpop
    for e in emmigrants[k]   # meta_pop[k] is the source, meta_pop[j] is the destination
      i = id[1]
      #println("e: ",e,"  i: ",i)
      #println("new emmigrant i: ",i,"  subpop_index:",k,"  num_attributes: ",num_attributes )
      vt[i] = deepcopy(vt[e])
      vt[i].fitness = fitness( vt[i].attributes, ideal )  
      #println("vt[",e,"]: ",vt[e])
      #println("vt[",i,"]: ",vt[i])
      Base.push!( new_emmigrants[j], i )
      id[1] += 1
    end
  end
  #println("new emmigrants: ",new_emmigrants)
  for j = 1:tr.num_subpops
    pop_after_deletion = Population[]
    #println("j: ",j,"  j%tr.num_subpops+1: ",j%tr.num_subpops+1,"  (j+tr.num_subpops-2)%tr.num_subpops+1: ",(j+tr.num_subpops-2)%tr.num_subpops+1)
    if neg_select  # use reverse proportional selection to delete elements by negative fitness
      pop_after_deletion = reverse_propsel(meta_pop[j],tr.ne,vt)
    else  # delete random elements to delete
      pop_after_deletion = meta_pop[j][1:(subpop_size-tr.ne)]
    end
    meta_pop[j] = append!( pop_after_deletion, new_emmigrants[j] )
  end
  #emmigrants  # perhaps should be the modified meta_pop
  meta_pop
end

function mutate_meta_pop!( meta_pop::PopList, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64}, mutation_stddev::Float64 )
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
      vt[meta_pop[j][i]] = mutate_variant( v_lists[j][i], mutation_stddev )
      vt[meta_pop[j][i]].fitness = fitness( vt[meta_pop[j][i]].attributes, ideal )
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
  Add a uniform random float between -move_range and +move_range to each ideal value.
"""
function move_optima( ideal::Vector{Float64}, move_range::Float64 )
  num_attributes = length(ideal)
  for k = 1:num_attributes
    move = 2.0*rand()*move_range - move_range
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
  for fm = mmeans
    if fm < opt_loss_cutoff
      count += 1
    end
  end
  count
end


#run_evolve(N,num_attributes,num_subpops, ngens,mutation_stddev,num_emmigrants,move_range,move_time_interval,opt_loss_cutoff,uniform_start)
function init()
  #include("types.jl")
  include("propsel.jl")
  include("fitness.jl")
  global N=64
  global num_subpops = 8
  global num_attributes = 4
  global ngens = 200
  global num_emmigrants = 2
  global mutation_stddev = 0.04
  global move_range = 0.1
  global move_time_interval = 5
  global opt_loss_cutoff = 0.2
  global horiz_select = false
  global uniform_start = false
end  
 
