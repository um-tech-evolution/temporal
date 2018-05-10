export mutate_meta_pop!, mutate_subpop!, mutate_attributes!

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
function mutate_meta_pop!( meta_pop::PopList, vt::Dict{Int64,variant_type}, ideal::Vector{Float64}, id::Vector{Int64}, tp::temporal_result_type  )
  num_subpops = length(meta_pop)
  subpop_size = length(meta_pop[1])
  #println("num_subpops: ",num_subpops,"  subpop_size: ",subpop_size)
  gen_innov_counts = generational_innovation_counts(0,0,0,0)
  for j = 1:num_subpops
    mutate_subpop!(meta_pop[j],vt,ideal,id,tp,gen_innov_counts)
  end
  return gen_innov_counts
end

function mutate_subpop!( subpop::Population, vt::Dict{Int64,variant_type}, ideal::Vector{Float64},
    id::Vector{Int64}, tp::temporal_result_type, innov_counts::generational_innovation_counts  )
  #println("starting mutate_subpop!  subpop: ",subpop)
  subpop_size = length(subpop)
  v_list = [vt[subpop[i]] for i = 1:subpop_size ]
  #innov_counts = generational_innovation_counts(0,0,0,0)
  fitnesses = zeros(Float64,subpop_size)
  for i = 1:subpop_size
    subpop[i] = id[1]
    id[1] += 1
    #vt[subpop[i]] = mutate_variant( v_list[i], tp[:mutStddev] )
    vt[subpop[i]] = deepcopy( v_list[i] )
    mutate_attributes!(  vt[subpop[i]].attributes, tp[:mutStddev], true )  # true means additive mutation instead of multiplictive
    fitnesses[i] = fitness( vt[subpop[i]].attributes, ideal, minFit=tp[:minFit], linfit_slope=tp[:linfit_slope] )
    vt[subpop[i]].fitness = fitnesses[i]
  end
  mmean = mean(fitnesses)
  #println("mutate_subpop: fitness mean: ",mmean)
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

#=
function mutate_variant( v::variant_type, mutStddev::Float64 )
  new_v = deepcopy(v)
  new_v.attributes = mutate_attributes( new_v.attributes, mutStddev, true )
  new_v
end
=#


function mutate_attributes!( attributes::Vector{Float64}, normal_stddev::Float64, additive_error::Bool )
  # Assumes that attributes can be modified in-place.
  #stddev = normal_stddev()   # Standard deviation of mutational perturbations
  #println("mutate attributes  normal_stddev: ",normal_stddev)
  #attributes = min(1.0,max(0.0,attributes+normal_stddev*randn(length(attributes))))
  #attributes = deepcopy(attributes)
  if additive_error
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
  else  # multiplicative copy error
    for i = 1:length(attributes)
      if attributes[i] <= 0.0
        #println("neg attribute: ",attributes[i])
        attributes[i] = 1.0e-6
      end
      multiplier = (1.0+normal_stddev*randn())
      #println("multiplier: ",multiplier)
      while multiplier <= 1.0e-6
        #println("neg multiplier")
        multiplier = (1.0+normal_stddev*randn())
      end
      attributes[i] *= multiplier
      if attributes[i] < 0.0
        #println("negative attribute with i=",i,": attribute: ",attributes[i])
      end
    end
  end
  #println("attributes: ",attributes)
  #return attributes
end

