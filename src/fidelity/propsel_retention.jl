# Model propsel_first in a population with 2 fitness levels and proportional selection.
# High fitness individuals have fitness 1+s, and low fitness individuals have fitness 1.
# Popsize is N.
# Mutation from high fitness to low fitness happens with probability 1-q.  
# There is no mutation in the other direction.
export count_high_low_fixed
using Distributions



#@everywhere function run_trial( N::Int64, q::Float64, s::Float64, ngens::Int64, init_high_freq::Float64, psel_first::Bool )
function run_trial( fd::fidelity_type )
  if fd.psel_first
    (g,k) = propsel_first( fd )
  else
    (g,k) = mutation_first( fd )
    #(g,k) = mutation_by_demotion( fd )
  end
  return (g,k)
end

function count_high_low_fixed( fd::fidelity_type )
  if fd.psel_first
    println("T: ",fd.T,"  N: ",fd.N,"  maxgens: ",fd.ngens,"  s: ",fd.s,"  q: ",fd.q,"  predicted p: ",(fd.s*fd.q+fd.q-1.0)/fd.s)
  else
    println("T: ",fd.T,"  N: ",fd.N,"  maxgens: ",fd.ngens,"  s: ",fd.s,"  q: ",fd.q,"  predicted p: ",(fd.s*fd.q+fd.q-1.0)/(fd.q*fd.s))
  end
  results = Tuple{Int64,Int64}[]
  #results = map( t->run_trial(deepcopy(fd)), [t for t in 1:T] )
  for t in 1:fd.T
    Base.push!( results, run_trial(deepcopy(fd) ) )
  end
  high_fixed = Dict{Int,Int}()   # maps number of generations to count of trials with high fixed in that number of generations
  low_fixed = Dict{Int,Int}()    # maps number of generations to count of trials with high fixed in that number of generations
  gen_reached = Dict{Int,Int}()  # maps k to the number of trials that reached the ngens limit with that k value
  #println("length(results): ",length(results))
  t=1
  for g_k in results
    g = g_k[1]
    k = g_k[2]
    #println("t: ",t,"  g: ",g,"  k: ",k)
    if k == fd.N
      if get(high_fixed,g,-1) == -1
        high_fixed[g] = 1
      else
        high_fixed[g] += 1
      end
    elseif k == 0
      if get(low_fixed,g,-1) == -1
        low_fixed[g] = 1
      else
        low_fixed[g] += 1
      end
    else
      if get(gen_reached,k,-1) == -1
        gen_reached[k] = 1
      else
        gen_reached[k] += 1
      end
    end
    t += 1
  end
  fd.high_fixed_count = sum_values(high_fixed)
  fd.low_fixed_count = sum_values(low_fixed)
  fd.gen_reached_count = sum_values(gen_reached)
  fd.average_high = average_k(gen_reached)
  #return (fd.high_fixed_count,fd.low_fixed_count,fd.gen_reached_count)
  return fd
end

@doc """ function propsel_first()
  Do up to ngens generations where a generation does proportional selection followed by mutation
"""
#@everywhere function propsel_first(fd.N::Int64, fd.q::Float64, s::Float64, ngens::Int64, init_high_frefd.q::Float64 )
function propsel_first( fd::fidelity_type )
  g = 0
  k = Int(round( fd.init_high_freq*fd.N ))  # k is the initial frequency of high fitness indivs
  done = false
  #println("g: ",g,"  k: ",k)
  while g <= fd.ngens && !done
    g += 1
    p = Float64(k)/Float64(fd.N)
    prob_high = (p+p*s)/(1.0+p*fd.s)
    h = rand(Binomial(fd.N,prob_high))  # h is the number of high fitness indivs after prop sel
    r = Float64(h)/fd.N                 # r is the proportion of high fitness indivs after prop sel
    k = rand(Binomial(fd.N,fd.q*r))        # k is the number of high fitness indivs after mutation
    #println("g: ",g,"  prob_high: ",prob_high,"  h: ",h,"  r: ",r,"  fd.q*r: ",fd.q*r,"  k: ",k)
    if k == 0 || k == fd.N
      done = true
    end
  end
  return (g,k)
end

@doc """ function mutation_first()
  Do up to ngens generations where a generation does proportional selection followed by mutation
"""
#@everywhere function mutation_first(fd.N::Int64, fd.q::Float64, s::Float64, ngens::Int64, init_high_frefd.q::Float64 )
function mutation_first( fd::fidelity_type )
  g = 0
  k = Int(round( fd.init_high_freq*fd.N ))  # k is the initial frequency of high fitness indivs
  done = false
  #println("g: ",g,"  k: ",k)
  while g <= fd.ngens && !done
    g += 1
    p = Float64(k)/Float64(fd.N)
    h = rand(Binomial(fd.N,fd.q*p))        # h is the number of high fitness indivs after mutation
    r = Float64(h)/Float64(fd.N)        # r is the proportion of high fitness indivs after mutation
    prob_high = (r+r*fd.s)/(1.0+r*fd.s)
    k = rand(Binomial(fd.N,prob_high))  # h is the number of high fitness indivs after prop sel
    #println("g: ",g,"  h: ",h,"  r: ",r,"  fd.q*p: ",fd.q*p,"  prob_high: ",prob_high,"  k: ",k)
    if k == 0 || k == fd.N
      done = true
    end
  end
  return (g,k)
end

@doc """  function mutation_by_demotion{)
  Implement mutation by deterministically demoting the fraction fd.q of high fitness individuals to low fitness.
"""
function mutation_by_demotion( fd::fidelity_type )
  g = 0
  k = Int(round( fd.init_high_freq*fd.N ))  # k is the initial frequency of high fitness indivs
  done = false
  #println("g: ",g,"  k: ",k)
  while g <= fd.ngens && !done
    g += 1
    p = Float64(k)/Float64(fd.N)
    prob_high = (p+p*s)/(1.0+p*fd.s)
    h = rand(Binomial(fd.N,prob_high))  # h is the number of high fitness indivs after prop sel
    k = Int(round(h*q))                 # p is the new frequency of high fitness individuals after deterministically demoting a fraction q
    println("g: ",g,"  prob_high: ",prob_high,"  h: ",h,"  p: ",p,"  k: ",k)
    if k == 0 || k == fd.N
      done = true
    end
  end
  return (g,k)
end

function average_k( ngens_reached::Dict{Int,Int} )
  if length( ngens_reached ) > 0
    sum_count = 0
    sum_k = 0
    for k in keys(ngens_reached)
      sum_count += ngens_reached[k]
      sum_k += k*ngens_reached[k]
    end
    return Float64(sum_k)/sum_count
  else
    return -1
  end
end

function sum_values( d::Dict{Int,Int} )
  sum_values = 0
  for k in keys(d)
    sum_values += d[k]
  end
  sum_values
end

