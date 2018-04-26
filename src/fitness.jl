export fitness
@doc """ euclidean_distance()
Euclidean distance from attributes to ideal
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

@doc """ sum_distance()
  Sum of absolute differences between ideal components and attribute components,
    except that if absolute difference is less than minFit, it is set to zero.
  count is the number of distance components that are greater than zero.
"""
function sum_distance( attributes::Vector{Float64}, ideal::Vector{Float64}, minFit::Float64 )
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  sum = 0.0
  count = 0
  for k = 1:length(attributes)
    a = abs(attributes[k]-ideal[k])
    a = a > 0.5 ? 1-a : a
    dis = 0.5-a
    if dis >= minFit
      # Don't modify dis
      count += 1
    else
      dis = 0.0
    end
    #println("k: ",k,"  attr: ",attributes[k],"  ideal: ",ideal[k],"  a: ",a,"  dis: ",dis)
    sum += dis
  end
  #println("fitness: attributes: ",attributes,"  ideal: ",ideal," fit: ",sum,"  count: ",count)
  #return sum, count
  return sum 
end

max_distance( attributes::Vector{Float64}, ideal::Vector{Float64} )= maximum(abs.(attributes-ideal))

@doc """ fitness()

"""
function fitness( attributes::Vector{Float64}, ideal::Vector{Float64}; 
      minFit::Float64=0.0, linfit_slope::Float64=1.0 )
  dis = euclidean_distance( attributes, ideal )
  fit = max( minFit, 0.5-linfit_slope*euclidean_distance( attributes, ideal ))
  #println("l fit: ",max( fit, minFit ))
  return max( fit, minFit )
end
