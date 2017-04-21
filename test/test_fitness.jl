export fitness, fit_counter
# Generate artificial fitness values for testing.

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

fit_counter = 0   # fit_counter global to the enclosing module
#=
function fit_init()
  global fit_first = true
  global fit_counter = 0
end
=# 

@doc """ fitness()
Artificially generated fitnesses that don't depend on the attributes and the ideal.
"""
function fitness( attributes::Vector{Float64}, ideal::Vector{Float64}; use_atan::Bool=false,
    min_fit::Float64=0.0, linear_fitness::Bool=false )
  global fit_counter
  fit_inc = 1.0
  if fit_counter % 2 == 0
    fit_counter += 1
    #println("test fit: fc: ",fit_counter,"  fit: ",min_fit)
    return min_fit
  else
    fit_counter += 1
    #println("test fit: fc: ",fit_counter,"  fit: ",min_fit+fit_inc)
    return min_fit+fit_inc
  end
end
  
  
