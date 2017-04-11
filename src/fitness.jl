export fitness
@doc """ dis()
Wraparoud Euclidean distance 
"""
function euclidean_distance( attributes::Vector{Float64}, ideal::Vector{Float64} )
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  sum = 0.0
  for k = 1:length(attributes)
    a = abs(attributes[k]-ideal[k])
    a = a > 0.5 ? 1-a : a
    sum += a^2
  end
  #println("fitness: attributes: ",attributes,"  ideal: ",ideal," fit: ",sqrt(sum))
  return sqrt(sum)
end

@doc """ fitness()
A strongly-peaked fitness function.
TODO:  add the parameters to the arguments.
"""
function fitness( attributes::Vector{Float64}, ideal::Vector{Float64}; use_atan::Bool=false,
    min_fit::Float64=0.0, linear_fitness::Bool=false )
  dis = euclidean_distance( attributes, ideal )
  if linear_fitness
    fit = max( min_fit, 0.5-euclidean_distance( attributes, ideal ))
    #println("l fit: ",fit)
    return fit
  end
  if use_atan
    fit = max(0.0, 1.0 - 2.0*atan(10.0*dis)/3.0)
    #println("a fit: ",fit)
  else
    inverse_scale = 20.0
    gpdf = Distributions.Gamma(1.0,1.0/inverse_scale)
    fit = pdf(gpdf,dis)/inverse_scale
    #println("g fit: ",fit)
  end
  return max( fit, min_fit )
end
