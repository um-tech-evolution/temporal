export fit_means_vars_stddevs_cfvars, attr_means_vars_stddevs_cfvars, fit_means 
@doc """ fit_means_vars_stddevs_cfvars()
means[i]  is the mean fitness of the individuals of subpops[i]
vars[i]   is the variance of the fitnesses of the individuals of subpops[i]
stds[i]   is the standard deviation of the fitnesses of the individuals of subpops[i] 
cfvars[i] is the coefficient of variation of the fitnesses of the individuals of subpops[i]
"""
function fit_means_vars_stddevs_cfvars( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  stds = [ std(map(fit,s)) for s in subpops ]
  cfvars = [ coef_var(map(fit,s)) for s in subpops ]
  return means, vars, stds, cfvars
end

function fit_means( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type} )
  fit(v) = variant_table[v].fitness
  fmeans = [ mean(map(fit,s)) for s in subpops ]
  return fmeans
end

function attr_means_vars_stddevs_cfvars( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type})
  num_attributes = length(variant_table[subpops[1][1]].attributes)
  (cumm_means, cumm_vars, cumm_stds, cumm_cfvars) = (0.0, 0.0, 0.0, 0.0)
  for j = 1:num_attributes
    #println("j: ",j,"  (i,s): ",[[(i,s) for i in s] for s in subpops])
    means = [ mean([variant_table[i].attributes[j] for i in s]) for s in subpops ]
    #println("means: ",means)
    vars = [ var([variant_table[i].attributes[j] for i in s]) for s in subpops ]
    #println("vars: ",vars)
    stds = [ std([variant_table[i].attributes[j] for i in s]) for s in subpops ]
    #println("stds: ",stds)
    cfvars = [ coef_var([variant_table[i].attributes[j] for i in s]) for s in subpops ]
    #println("cfvars: ",cfvars)
    (cumm_means, cumm_vars, cumm_stds, cumm_cfvars) = 
        (cumm_means, cumm_vars, cumm_stds, cumm_cfvars) .+ (means, vars, stds, cfvars )
  end
  return (cumm_means, cumm_vars, cumm_stds, cumm_cfvars)./num_attributes
end


@doc """ attr_vars()
Returns ave_vars, where ave_vars[s] is the mean (over attributes) of the variance of each attribute, where the
  variance is taken over the members of subpopulation s.
"""
function attr_vars( subpops::PopList, variant_table::Dict{Int64,temporal_variant_type} )
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

function coef_var( lst )
  return std(lst)/mean(lst)
end

