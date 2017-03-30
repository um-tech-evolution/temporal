export fitness
@doc """ dis()
Euclidean distance
"""
function euclidean_distance( attributes::Vector{Float64}, ideal::Vector{Float64} )
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  sum = 0.0
  for k = 1:length(attributes)
    sum += ( attributes[k] - ideal[k] )^2
  end
  #println("fitness: attributes: ",attributes,"  ideal: ",ideal," fit: ",sqrt(sum))
  return sqrt(sum)
end

@doc """ fitness()
A strongly-peaked fitness function.
TODO:  add the parameters to the arguments.
"""
function fitness( attributes::Vector{Float64}, ideal::Vector{Float64}; use_atan::Bool=false )
  dis = euclidean_distance( attributes, ideal )
  if use_atan
    return max(0.0, 1.0 - 2.0*atan(10.0*dis)/3.0)
  else
    inverse_scale = 20.0
    gpdf = Distributions.Gamma(1.0,1.0/inverse_scale)
    return pdf(gpdf,dis)/inverse_scale
  end
end
